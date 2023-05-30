# Figure S10: Copy number filtering by PTATO
library(ggplot2)
library(tidyverse)
library(devtools)
library(magrittr)
library(MutationalPatterns)
library(reshape2)

input_dir <- "/path/to/MendeleyData/directory/"
output_dir <- "/path/to/output/directory/"
if(dir.exists(output_dir) == F) { dir.create(output_dir)}

# S10A = Recurrency PTA readcounts
readCount_files <- list("IBFM35-DX2BM-AMLBULK" = paste(input_dir, "IBFM35/PTATO/intermediate/cnvs/cobalt/IBFM35/IBFM35-DX2BM-MSCBULK/IBFM35-DX2BM-AMLBULK.cobalt.ratio.tsv", sep = ""),
                        "IBFM35-DX2BM-HSCPTAP1D9" = paste(input_dir, "IBFM35/PTATO/intermediate/cnvs/cobalt/IBFM35/IBFM35-DX2BM-MSCBULK/IBFM35-DX2BM-HSCPTAP1D9.cobalt.ratio.tsv", sep = ""),
                        "IBFM35-DX2BM-HSCPTAP1E9" = paste(input_dir, "IBFM35/PTATO/intermediate/cnvs/cobalt/IBFM35/IBFM35-DX2BM-MSCBULK/IBFM35-DX2BM-HSCPTAP1E9.cobalt.ratio.tsv", sep = ""),
                        "IBFM35-DX2BM-HSCPTAP1G9" = paste(input_dir, "IBFM35/PTATO/intermediate/cnvs/cobalt/IBFM35/IBFM35-DX2BM-MSCBULK/IBFM35-DX2BM-HSCPTAP1G9.cobalt.ratio.tsv",sep = ""))

readCounts <- data.frame()
readCounts_summary <- data.frame()
for(sample in names(readCount_files)){
  print(sample)
  print(readCount_files[[sample]])
  # Read the first 10000 lines / sample
  readCounts_sample <- read.delim(readCount_files[[sample]])
  readCounts_sample$Sample <- sample
  readCounts_summary <- rbind(readCounts_summary, data.frame(Sample = sample, Total_readcounts = sum(readCounts_sample$tumorReadCount)))
  readCounts <- rbind(readCounts, readCounts_sample)
}

readCounts_summary$Normalization <- readCounts_summary$Total_readcounts / mean(readCounts_summary$Total_readcounts)

#readCounts$CopyNumber <- readCounts$

head(readCounts)

# Select a 50kb region on chr1 as an example
readCounts_plot <- readCounts[readCounts$chromosome == 1 & readCounts$position > 7868251 & readCounts$position < (7868251+5e4),]
readCounts_plot <- merge(readCounts_plot, readCounts_summary, by = "Sample", all.x = T)
readCounts_plot$ReadCounts_normalized <- readCounts_plot$tumorReadCount/readCounts_plot$Normalization
readCounts_plot <- readCounts_plot[,c("Sample", "chromosome", "position", "ReadCounts_normalized")]
readCounts_plot$Sample[readCounts_plot$Sample == "IBFM35-DX2BM-HSCPTAP1D9"] <- paste(readCounts_plot$Sample[readCounts_plot$Sample == "IBFM35-DX2BM-HSCPTAP1D9"] , "\n(IBFM35_1)")
readCounts_plot$Sample[readCounts_plot$Sample == "IBFM35-DX2BM-HSCPTAP1E9"] <- paste(readCounts_plot$Sample[readCounts_plot$Sample == "IBFM35-DX2BM-HSCPTAP1E9"] , "\n(IBFM35_2)")
readCounts_plot$Sample[readCounts_plot$Sample == "IBFM35-DX2BM-HSCPTAP1G9"] <- paste(readCounts_plot$Sample[readCounts_plot$Sample == "IBFM35-DX2BM-HSCPTAP1G9"] , "\n(IBFM35_3)")

ggplot(readCounts_plot, aes(x = position, y= ReadCounts_normalized, col = Sample)) + geom_line(linewidth = 0.4, alpha = 0.8) +
  #scale_color_brewer(palette = "Dark2") +
  scale_color_manual(values = c( "darkorange","#377EB8","#4DAF4A","#984EA3")) +
  labs(y = "Readcounts per 1kb bin", x = "Genomic position (chr1)") +
  theme_classic(base_size = 6) + scale_x_continuous(expand = c(0,0)) 

ggsave(paste(output_dir, "FigS10A_IBFM35_CoverageChr1.pdf"), width = 13, height = 5, units = "cm")

# S10B = Cosine similarities in (COBALT) coverage patterns PTA samples
trl_tbl = tibble("sample_name" = c("IBFM35-DX2BM-AMLBULK", "IBFM35-DX2BM-HSCPTAP1D9", "IBFM35-DX2BM-HSCPTAP1E9", "IBFM35-DX2BM-HSCPTAP1G9"), 
                 "label" = c("AML Bulk", "IBFM35_1", "IBFM35_2", "IBFM35_3"))

# Find samples
cobalt_fnames = list.files( paste(input_dir, "IBFM35/PTATO/intermediate/cnvs/cobalt/IBFM35/IBFM35-DX2BM-MSCBULK/", sep = ""), pattern = ".ratio.tsv", full.names = T)
cobalt_fnames = str_subset(cobalt_fnames, "HSCPTAP1D9|HSCPTAP1E9|HSCPTAP1G9|AMLBULK")
sample_names = cobalt_fnames %>% 
  basename() %>% 
  str_remove(".cobalt.ratio.tsv")

# Read data
readcount_l = purrr::map(cobalt_fnames, ~read_tsv(.x, col_select = c(chromosome, position, tumorReadCount))) %>% 
  purrr::map2(sample_names, ~set_colnames(.x, c("chromosome", "position", .y)))

# Combine all samples in a single tibble
read_count_tbl = purrr::reduce(readcount_l, inner_join, by = c("chromosome", "position"))

# Calculate cosine similarity
cos_sim_m = read_count_tbl %>% 
  dplyr::select(-chromosome, -position) %>% 
  cos_sim_matrix(., .)

# Replace sample names with labels for manuscript
labels = tibble("sample_name" = colnames(cos_sim_m)) %>% 
  left_join(trl_tbl, by = "sample_name") %>% 
  pull(label)
colnames(cos_sim_m) = labels
rownames(cos_sim_m) = labels

# Calculate mean/median cosine similarities
cos_sim_m_ptaonly = cos_sim_m[2:4, 2:4]
cossims_pta_only = cos_sim_m_ptaonly[lower.tri(cos_sim_m_ptaonly)]
mean(cossims_pta_only)
quantile(cossims_pta_only)

cossims_with_bulk = cos_sim_m[1, 2:4]
mean(cossims_with_bulk)
quantile(cossims_with_bulk)

# function
plot_cosine_heatmap <- function(cos_sim_matrix, col_order = NA, row_order = NA, cluster_rows = TRUE,
                                cluster_cols = FALSE, method = "complete", plot_values = FALSE) {
  # check explained argument
  if (!inherits(cos_sim_matrix, "matrix")) {
    stop("cos_sim_matrix must be a matrix")
  }
  # matrix should have row and colnames
  if (length(colnames(cos_sim_matrix)) == 0) {
    stop("cos_sim_matrix is missing colnames")
  }
  if (length(rownames(cos_sim_matrix)) == 0) {
    stop("cos_sim_matrix is missing rownames")
  }
  
  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  Cosine.sim <- Signature <- Sample <- x <- y <- xend <- yend <- NULL
  
  # If cluster_rows is TRUE perform clustering. Else use supplied row_order or
  # the current column order.
  if (!is.na(row_order) & cluster_rows == TRUE) {
    stop("row_order can only be provided when cluster_rows is FALSE", call. = FALSE)
  } else if (!is.na(row_order)) {
    # check row_order argument
    if (!inherits(row_order, "character")) {
      stop("row_order must be a character vector", call. = FALSE)
    }
    if (length(row_order) != nrow(cos_sim_matrix)) {
      stop("row_order must have the same length as the number of
          samples in the explained matrix", call. = FALSE)
    }
  } else if (cluster_rows == TRUE) {
    # cluster samples based on euclidean distance between relative contribution
    hc.sample <- hclust(dist(cos_sim_matrix), method = method)
    # order samples according to clustering
    row_order <- rownames(cos_sim_matrix)[hc.sample$order]
    
    dhc <- as.dendrogram(hc.sample)
    # rectangular lines
    ddata <- ggdendro::dendro_data(dhc, type = "rectangle")
    # plot dendrogram of hierachical clustering
    dendrogram_rows <- ggplot(ggdendro::segment(ddata)) +
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
      coord_flip() +
      scale_y_reverse(expand = c(0.2, 0)) +
      ggdendro::theme_dendro()
  }
  else {
    row_order <- rownames(cos_sim_matrix)
  }
  
  
  # If cluster_cols is TRUE perform clustering. Else use supplied col_order or
  # the current column order.
  if (!is.na(col_order) & cluster_cols == TRUE) {
    stop("col_order can only be provided when cluster_cols is FALSE", call. = FALSE)
  } else if (!is.na(col_order)) {
    # check col_order argument
    if (!inherits(col_order, "character")) {
      stop("col_order must be a character vector", call. = FALSE)
    }
    if (length(col_order) != ncol(cos_sim_matrix)) {
      stop("col_order must have the same length as the number of 
          signatures in the explained matrix", call. = FALSE)
    }
  } else if (cluster_cols == TRUE) {
    # Cluster cols
    hc.sample2 <- cos_sim_matrix %>%
      t() %>%
      dist() %>%
      hclust(method = method)
    col_order <- colnames(cos_sim_matrix)[hc.sample2$order]
    
    dhc <- as.dendrogram(hc.sample2)
    # rectangular lines
    ddata <- ggdendro::dendro_data(dhc, type = "rectangle")
    # plot dendrogram of hierachical clustering
    dendrogram_cols <- ggplot(ggdendro::segment(ddata)) +
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
      ggdendro::theme_dendro() +
      scale_y_continuous(expand = c(0.2, 0))
  } else {
    col_order <- colnames(cos_sim_matrix)
  }
  
  
  # Make matrix long and set factor levels, to get the correct order for plotting.
  cos_sim_matrix.m <- cos_sim_matrix %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Sample") %>%
    tidyr::pivot_longer(-Sample, names_to = "Signature", values_to = "Cosine.sim") %>%
    dplyr::mutate(
      Signature = factor(Signature, levels = col_order),
      Sample = factor(Sample, levels = row_order)
    )
  
  # plot heatmap
  heatmap <- ggplot(cos_sim_matrix.m, aes(x = Signature, y = Sample, fill = Cosine.sim, order = Sample)) +
    geom_raster() +
    scale_fill_distiller(palette = "YlGnBu", direction = 1, name = "Cosine \nsimilarity", limits = c(0, 1.000000001)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = NULL, y = NULL)
  
  # if plot_values is TRUE, add values to heatmap
  if (plot_values == TRUE) {
    heatmap <- heatmap + geom_text(aes(label = round(Cosine.sim, 2)), size = 1.6 , col = "white")
  }
  
  # Add dendrogram depending on the clustering of the rows and the columns.
  if (cluster_rows == TRUE & cluster_cols == TRUE) {
    empty_fig <- ggplot() +
      theme_void()
    plot_final <- cowplot::plot_grid(empty_fig, dendrogram_cols, dendrogram_rows, heatmap,
                                     align = "hv", axis = "tblr", rel_widths = c(0.3, 1), rel_heights = c(0.3, 1)
    )
  }
  else if (cluster_rows == TRUE & cluster_cols == FALSE) {
    # combine plots
    plot_final <- cowplot::plot_grid(dendrogram_rows, heatmap, align = "h", rel_widths = c(0.3, 1))
  } else if (cluster_rows == FALSE & cluster_cols == TRUE) {
    plot_final <- cowplot::plot_grid(dendrogram_cols, heatmap, align = "v", rel_heights = c(0.3, 1)) +
      # reverse order of the samples such that first is up
      ylim(rev(levels(factor(cos_sim_matrix.m$Sample))))
  } else {
    plot_final <- heatmap +
      # reverse order of the samples such that first is up
      ylim(rev(levels(factor(cos_sim_matrix.m$Sample))))
  }
  
  return(plot_final)
}

# Plot cosine similarity
cossim_fig = plot_cosine_heatmap(cos_sim_m, cluster_rows = FALSE, cluster_cols = FALSE, plot_values = TRUE) + 
  theme_classic(base_size = 6) + theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1),legend.key.height= unit(0.3, 'cm'),
                                       legend.key.width= unit(0.4, 'cm'))
ggsave(paste(output_dir, "S10B_IBFM35_cobalt_cossims.pdf", sep = ""), cossim_fig, width = 4.5, height = 5, units = "cm")



### Figure S10C-S10H

cobalt_filtered_file <- paste(input_dir, "IBFM35/PTATO/intermediate/cnvs/readCounts/IBFM35-DX2BM-HSCPTAP1G9/IBFM35-DX2BM-HSCPTAP1G9.readcounts.filtered.1kb.txt", sep = "")
cobalt_filtered <- read.delim(cobalt_filtered_file)

cobalt_100kb_file <- paste(input_dir, "IBFM35/PTATO/intermediate/cnvs/readCounts/IBFM35-DX2BM-HSCPTAP1G9/IBFM35-DX2BM-HSCPTAP1G9.readcounts.filtered.100kb.txt", sep = "")
cobalt_100kb <- read.delim(cobalt_100kb_file)

segments_file <- paste(input_dir, "IBFM35/PTATO/intermediate/cnvs/readCounts/IBFM35-DX2BM-HSCPTAP1G9/IBFM35-DX2BM-HSCPTAP1G9.readcounts.segments.txt", sep = "")
segments <- read.delim(segments_file)

copyNumberPlot <- function(input_data, column_y, y_label = "Copy number"){
  
  plot <- 
    ggplot(input_data, aes(x = position, y = input_data[,column_y])) + 
    geom_point(size = 0.02, alpha = 0.2) + 
    facet_grid(.~chromosome, space = "free_x", scales = "free_x",switch = "x") +
    theme_classic(base_size = 5) +
    coord_cartesian(ylim = c(0,6)) +
    labs(y = y_label) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.minor=element_blank(),
          panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_line(colour="grey90", linetype=2, linewidth=0.2),
          panel.spacing = unit(0.2, "mm"),
          panel.background=element_blank(),
          panel.border=element_blank(),
          strip.background = element_blank(),
          strip.placement = "outside")
  
  return(plot)
}

raw_plot <- copyNumberPlot(input_data = cobalt_filtered[cobalt_filtered$chromosome %in% c(1:3),],
                              column_y = "copyNumber_raw", y_label = "Copy number / 1kb\n(Raw)")
raw_plot
ggsave(paste(output_dir, "S10D_IBFM35-DX2BM-HSCPTAP1G9_CNraw.png", sep = ""), dpi = 600, width = 6, height = 2, units = "cm")

smooth_plot <- copyNumberPlot(input_data = cobalt_filtered[cobalt_filtered$chromosome %in% c(1:3) & cobalt_filtered$FILTER == "PASS",],
                              column_y = "copyNumber_smooth", y_label = "Copy number / 1kb\n(Smoothened)")
smooth_plot
ggsave(paste(output_dir, "S10E_IBFM35-DX2BM-HSCPTAP1G9_CNsmooth.png", sep = ""), dpi = 600, width = 6, height = 2, units = "cm")

ggplot(cobalt_100kb[cobalt_100kb$Chromosome %in% c(1:3),], aes(x = (start.pos+end.pos)/2, y = medianReadCount)) + 
  geom_point(size = 0.02, alpha = 0.2) + 
  facet_grid(.~Chromosome, space = "free_x", scales = "free_x",switch = "x") +
  theme_classic(base_size = 5) +
  coord_cartesian(ylim = c(0,6)) +
  labs(y = "Copy number\n/ 100kb") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="grey90", linetype=2, linewidth=0.2),
        panel.spacing = unit(0.2, "mm"),
        panel.background=element_blank(),
        panel.border=element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside")
ggsave(paste(output_dir, "S10F_IBFM35-DX2BM-HSCPTAP1G9_CNbinned.png", sep = ""), dpi = 600, width = 6, height = 2, units = "cm")

colnames(segments)[1] <- "Chromosome"
segments$Segment <- paste("#", 1:nrow(segments), sep = "")
ggplot(cobalt_100kb[cobalt_100kb$Chromosome %in% c(1:3),], aes(x = (start.pos+end.pos)/2, y = medianReadCount)) + 
  geom_point(size = 0.02, alpha = 0.2) + 
  geom_segment(data = segments[segments$Chromosome %in% c(1:3),], aes(x = start.pos, xend = end.pos, y = mean, yend = mean), col = "red") +
  geom_text(data = segments[segments$Chromosome %in% c(1:3),], aes(x = (start.pos + end.pos) / 2, label = Segment, y = 5), size = 1.5) +
  facet_grid(.~Chromosome, space = "free_x", scales = "free_x",switch = "x") +
  theme_classic(base_size = 5) +
  coord_cartesian(ylim = c(0,6)) +
  labs(y = "Copy number\n/ 100kb") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="grey90", linetype=2, linewidth=0.2),
        panel.spacing = unit(0.2, "mm"),
        panel.background=element_blank(),
        panel.border=element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside")
ggsave(paste(output_dir, "S10G_IBFM35-DX2BM-HSCPTAP1G9_CNsegmented.png", sep = ""), dpi = 600, width = 6, height = 2, units = "cm")


cobalt_filtered2 <- cobalt_filtered
cobalt_filtered2$copyNumber_deviation <- abs(2-cobalt_filtered2$copyNumber_smooth)

cobalt_normal_segments <- cobalt_filtered2[cobalt_filtered2$FILTER == "PASS",]
cobalt_normal_segments2 <- cobalt_normal_segments[cobalt_normal_segments$copyNumber_deviation < as.numeric(quantile(cobalt_normal_segments$copyNumber_deviation, 0.25, na.rm = T)),]
cobalt_normal_segments2$segment <- "Normal"


cobalt_filtered2$segment <- paste("#", cobalt_filtered2$segment, sep ="")
cobalt_filtered3 <- rbind(cobalt_filtered2, cobalt_normal_segments2)


cobalt_filtered3_m <- melt(cobalt_filtered3[ ,c("copyNumber_raw", "ponCopyNumber_mean", "segment")], id.vars = "segment")
cobalt_filtered3_m$segment <- factor(cobalt_filtered3_m$segment)

cobalt_filtered3_m$Label <- ifelse(cobalt_filtered3_m$variable == "copyNumber_raw", "Sample", "PON")

quantile(cobalt_normal_segments2$copyNumber_raw)

ggplot(cobalt_filtered3_m[cobalt_filtered3_m$segment %in% c(paste("#",1:9, sep =""), "Normal"),], aes(x = Label, y = value, fill = Label)) + 
  annotate("rect", ymin = as.numeric(quantile(cobalt_normal_segments2$copyNumber_raw,0.25)),
           ymax =  as.numeric(quantile(cobalt_normal_segments2$copyNumber_raw,0.75)), xmin = -Inf, xmax = Inf,
           alpha = .1,fill = "blue") +
  
  geom_boxplot(outlier.shape = NA, linewidth = 0.2) + facet_grid(.~segment) + coord_cartesian(ylim = c(0,6)) +
  scale_fill_manual(values = c( "#D55E00", "#009E73")) +
  ggtitle("Copy Number Segments") +
  
  labs(y = "Copy number\n/ 1kb", fill = "Read counts") +
  theme_bw(base_size = 5) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="grey90", linetype=2, linewidth=0.2),
        panel.spacing = unit(0.2, "mm"),
        legend.key.width = unit(0.25, "cm"), plot.title = element_text(hjust = 0.5, size = 5))
ggsave(paste(output_dir, "S10H_IBFM35-DX2BM-HSCPTAP1G9_segments.pdf", sep = ""), dpi = 300, width = 7.5, height = 2.2, units = "cm")

### 

cobalt_filtered_list <- list("IBFM35-DX2BM-HSCPTAP1G9" = cobalt_filtered)
cobalt_filtered_list[["IBFM35-DX2BM-HSCPTAP1D9"]] <- read.delim( paste(input_dir, "IBFM35/PTATO/intermediate/cnvs/readCounts/IBFM35-DX2BM-HSCPTAP1D9/IBFM35-DX2BM-HSCPTAP1D9.readcounts.filtered.1kb.txt", sep = ""))
cobalt_filtered_list[["IBFM35-DX2BM-HSCPTAP1E9"]] <- read.delim( paste(input_dir, "IBFM35/PTATO/intermediate/cnvs/readCounts/IBFM35-DX2BM-HSCPTAP1E9/IBFM35-DX2BM-HSCPTAP1E9.readcounts.filtered.1kb.txt", sep = ""))

cobalt_filtered_list <- lapply(cobalt_filtered_list, function(x) x[which(as.vector(x$chromosome) %in%  c(1:22)),])

cobalt_100kb_list <- list("IBFM35-DX2BM-HSCPTAP1G9" = cobalt_100kb)
cobalt_100kb_list[["IBFM35-DX2BM-HSCPTAP1D9"]] <- read.delim( paste(input_dir, "IBFM35/PTATO/intermediate/cnvs/readCounts/IBFM35-DX2BM-HSCPTAP1D9/IBFM35-DX2BM-HSCPTAP1D9.readcounts.filtered.100kb.txt", sep = ""))
cobalt_100kb_list[["IBFM35-DX2BM-HSCPTAP1E9"]] <- read.delim( paste(input_dir, "IBFM35/PTATO/intermediate/cnvs/readCounts/IBFM35-DX2BM-HSCPTAP1E9/IBFM35-DX2BM-HSCPTAP1E9.readcounts.filtered.100kb.txt", sep = ""))

cobalt_100kb_list <- lapply(cobalt_100kb_list, function(x) x[which(as.vector(x$Chromosome) %in%  c(1:22)),])

segments_list <-  list("IBFM35-DX2BM-HSCPTAP1G9" = read.delim( paste(input_dir, "IBFM35/PTATO/intermediate/cnvs/readCounts/IBFM35-DX2BM-HSCPTAP1G9/IBFM35-DX2BM-HSCPTAP1G9.readcounts.segments.txt", sep = "")))
segments_list[["IBFM35-DX2BM-HSCPTAP1D9"]] <- read.delim( paste(input_dir, "IBFM35/PTATO/intermediate/cnvs/readCounts/IBFM35-DX2BM-HSCPTAP1D9/IBFM35-DX2BM-HSCPTAP1D9.readcounts.segments.txt", sep = ""))
segments_list[["IBFM35-DX2BM-HSCPTAP1E9"]] <- read.delim( paste(input_dir, "IBFM35/PTATO/intermediate/cnvs/readCounts/IBFM35-DX2BM-HSCPTAP1E9/IBFM35-DX2BM-HSCPTAP1E9.readcounts.segments.txt", sep = ""))
segments_list <- lapply(segments_list, function(x) x[which(as.vector(x$chrom) %in%  c(1:22)),])

cobalt_overview <- data.frame()
for(sample in names(cobalt_filtered_list)){
  print(sample)
  cobalt_overview_sample <- rbind(data.frame(Sample = sample, Variable = "Raw",  value =  cobalt_filtered_list[[sample]]$copyNumber_raw),
                           data.frame(Sample = sample, Variable = "Smooth",  value = cobalt_filtered_list[[sample]]$copyNumber_smooth),
                           data.frame(Sample = sample, Variable = "Binned",  value = cobalt_100kb_list[[sample]]$medianReadCount),
                           data.frame(Sample = sample, Variable = "Segments",  value = segments_list[[sample]]$mean))
  cobalt_overview <- rbind(cobalt_overview, cobalt_overview_sample)
}

cobalt_overview$Variable <- factor(cobalt_overview$Variable, levels = c("Raw","Smooth","Binned","Segments"))
cobalt_overview$Sample2 <- "IBFM35_1"
cobalt_overview$Sample2[cobalt_overview$Sample == "IBFM35-DX2BM-HSCPTAP1E9"] <- "IBFM35_2"
cobalt_overview$Sample2[cobalt_overview$Sample == "IBFM35-DX2BM-HSCPTAP1G9"] <- "IBFM35_3"

ggplot(cobalt_overview, aes(  y = value, fill = Variable, x  = Variable)) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.2) + 
  facet_grid(.~Sample2, space = "free", scales = "free") +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(ylim = c(0,5)) +
  theme_classic(base_size = 6) +
  theme(panel.border =  element_rect(color = "black", fill = NA), strip.background = element_rect(fill = "gray90"), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(y = "Copy number per window/segment", fill = "Windows", x = "Copy number windows")

ggsave(paste(output_dir, "S10I_CopyNumber_Boxplots.pdf", sep = ""), width = 9, height = 4, units = "cm")




cobalt_sd_overview <- data.frame()
for(sample in names(cobalt_filtered_list)){
  print(sample)
  cobalt_sd_overview_sample <- rbind(data.frame(Sample = sample, Variable = "Raw",  value =  sd(cobalt_filtered_list[[sample]]$copyNumber_raw, na.rm = T)),
                                  data.frame(Sample = sample, Variable = "Smooth",  value = sd(cobalt_filtered_list[[sample]]$copyNumber_smooth[!is.infinite(cobalt_filtered_list[[sample]]$copyNumber_smooth)], na.rm = T)),
                                  data.frame(Sample = sample, Variable = "Binned",  value = sd(cobalt_100kb_list[[sample]]$medianReadCount, na.rm = T)),
                                  data.frame(Sample = sample, Variable = "Segments",  value = sd(segments_list[[sample]]$mean, na.rm = T)))
  cobalt_sd_overview <- rbind(cobalt_sd_overview, cobalt_sd_overview_sample)
}

cobalt_sd_overview$Variable <- factor(cobalt_sd_overview$Variable, levels = c("Raw","Smooth","Binned","Segments"))
cobalt_sd_overview$Sample2 <- "IBFM35_1"
cobalt_sd_overview$Sample2[cobalt_sd_overview$Sample == "IBFM35-DX2BM-HSCPTAP1E9"] <- "IBFM35_2"
cobalt_sd_overview$Sample2[cobalt_sd_overview$Sample == "IBFM35-DX2BM-HSCPTAP1G9"] <- "IBFM35_3"

ggplot(cobalt_sd_overview, aes( x = Variable, y = value, fill = Variable)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  facet_grid(.~Sample2, space = "free", scales = "free") +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(ylim = c(0,3)) +
  theme_classic(base_size = 6) +
  theme(panel.border =  element_rect(color = "black", fill = NA), strip.background = element_rect(fill = "gray90"), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.key.size = unit(0.4, "cm")) +
  labs(y = "Copy number\nStandard deviation", fill = "Windows", x = "Copy number windows")

ggsave(paste(output_dir, "S10J_CopyNumber_Standarddeviations.pdf", sep = ""), width = 9, height = 4, units = "cm")



