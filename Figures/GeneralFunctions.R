# This script contains several general functions that are used in multiple scripts

library(dplyr)
library(ggplot2)

# Color settings from MutationalPatterns
COLORS7 <- c(
  "#2EBAED", "#000000",
  "#E98C7B", "#DE1C14","#D4D2D2", "#ADCC54",
  "#F0D0CE"
)

# Color settings from MutationalPatterns
COLORS6 <- c(
  "#2EBAED", "#000000",
  "#DE1C14","#D4D2D2", "#ADCC54",
  "#F0D0CE"
)


# This function selects the PTAprobCutoffs from the header of a FILTERED PTATO VCF
getPTAprobCutoff <- function(VCF_file){
  VCF_header <- scanVcfHeader(VCF_file)
  PTAprobCutoff_raw <- meta(VCF_header)$PTAprobCutoff[1,1]
  PTAprobCutoff_A  <- gsub(pattern = "\"", replacement = "", x = PTAprobCutoff_raw, fixed=TRUE)
  PTAprobCutoff_B <- gsub(pattern = "\\(|\\)|- ", replacement = "", x = PTAprobCutoff_A)
  PTAprobCutoffs <- unlist(strsplit(PTAprobCutoff_B, split = " "))
  names(PTAprobCutoffs) <- c("PTAprobCutoff", "Walker", "Cossim")
  return(PTAprobCutoffs)
}


# This function can be used to transform a VariantAnnotation vcf to a MutationalPatterns grangeslist (human only)
VCF_to_GR <- function(VCF, chromosomes = c(1:22, "X", "Y"), add_chr = TRUE){
  GR <- rowRanges(VCF)
  GR <- GR[which(seqnames(GR) %in% chromosomes),]
  seqlevels(GR) <- seqlevels(GR)[1:24]
  # Mutationalpatterns usually uses chr
  if(add_chr == TRUE){
    seqlevels(GR) <- paste("chr", seqlevels(GR), sep = "")
  } 
  return(GR)
}

# Variants with a VAF<VAF_threshold will be flagged with FAIL_VAF in the FILTER field.
# Variants with PTAprob>PTAprobCutoff will be flagged with FAIL in the FILTER field.
readPTATOvcf <- function(vcf_unfiltered, 
                         vcf_filtered = "",
                         cutoff = 1, # 1 = mean, 2 = Lira, 3 = Cossim
                         type = "snv",
                         VAF_threshold = 0.25,
                         chromosomes = c(1:22, "X", "Y"),
                         genome = "hg38"){
  
  vcf <- readVcf(vcf_unfiltered, genome = genome)
  vcf <- vcf[as.vector(seqnames(vcf)) %in% chromosomes,]
  GR <- VCF_to_GR(vcf, chromosomes = chromosomes)
  
  if(type != "indel"){
    PTAprobCutoff <- getPTAprobCutoff(vcf_filtered)
  } else {
    # indels that are removed have a PTAprob of 1
    PTAprobCutoff <- 0.5
  }
  
  GR$VAF <- as.numeric(geno(vcf)$VAF)
  GR$PTAprob <- as.numeric(as.data.frame(geno(vcf)$PTAprob)[,1])
  GR$PTAprobCutoff <- PTAprobCutoff[1]
  GR$PTAprobCutoff_Walker <- PTAprobCutoff[2]
  GR$PTAprobCutoff_Cossim <- PTAprobCutoff[3]
  
  GR$FILTER[GR$VAF < VAF_threshold] <- "FAIL_VAF"
  GR$FILTER[GR$PTAprob > PTAprobCutoff[cutoff]] <- "FAIL"
  
  return(GR)
}

# This function can be used to adjust the background colors of the facet strips
# plot is a ggplot object. use pdf() to save the plot, ggsave does not work
plot_adjusted_strips <- function(plot, strip_colors){
  plot_table <- ggplot_gtable(ggplot_build(plot))
  strip_both <- which(grepl('strip-', plot_table$layout$name))
  fills <- strip_colors
  k <- 1
  for (i in strip_both) {
    j <- which(grepl('rect', plot_table$grobs[[i]]$grobs[[1]]$childrenOrder))
    plot_table$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  plot(plot_table)
}


### 96-profiles
plot_96_profile_small <- function(mut_matrix, colors = NA, ymax = 0.2, condensed = FALSE, add_total= T) {
  
  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  freq <- full_context <- substitution <- context <- NULL
  
  # Check color vector length
  # Colors for plotting
  if (is.na(colors)) {
    colors <- c(
      "#2EBAED", "#000000", "#DE1C14",
      "#D4D2D2", "#ADCC54", "#F0D0CE"
    )
  }
  if (length(colors) != 6) {
    stop("Provide colors vector with length 6", call. = FALSE)
  }
  
  total_mutations <- data.frame(sample = colnames(mut_matrix), total = colSums(mut_matrix))
  total_mutations$total <- round(total_mutations$total, 0)
  
  # Make contribution relative
  norm_mut_matrix <- apply(mut_matrix, 2, function(x) x / sum(x))
  
  if(add_total == TRUE){
    colnames(norm_mut_matrix) <- paste(colnames(norm_mut_matrix),"\n(n=", total_mutations$total, ")", sep = "")
  }
  
  # Get substitution and context from rownames and make long.
  tb <- norm_mut_matrix %>%
    as.data.frame() %>%
    tibble::rownames_to_column("full_context") %>%
    dplyr::mutate(
      substitution = stringr::str_replace(full_context, "\\w\\[(.*)\\]\\w", "\\1"),
      context = stringr::str_replace(full_context, "\\[.*\\]", "\\.")
    ) %>%
    dplyr::select(-full_context) %>%
    tidyr::pivot_longer(c(-substitution, -context), names_to = "sample", values_to = "freq") %>% 
    dplyr::mutate(sample = factor(sample, levels = unique(sample)))
  
  
  
  # Change plotting parameters based on whether plot should be condensed.
  if (condensed == TRUE) {
    width <- 1
    spacing <- 0
  } else {
    width <- 0.6
    spacing <- 0.5
  }
  
  # Create figure
  plot <- ggplot(data = tb, aes(x = context, y = freq, 
                                fill = substitution, width = width)) +
    geom_bar(stat = "identity", colour = "white", size = .1) +
    geom_hline(yintercept = 0, size = 0.5) +
    scale_fill_manual(values = colors) +
    facet_grid(sample ~ substitution) +
    ylab("Relative contribution") +
    coord_cartesian(ylim = c(0, ymax)) +
    scale_y_continuous(breaks = seq(0, ymax, ymax/2), expand = c(0,0)) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    guides(fill = "none") +
    theme_classic(base_size = 6) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(vjust = 1),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 4),
      panel.grid.major.x = element_blank(),
      panel.spacing.y = unit(0.5, "lines"),
      panel.spacing.x = unit(spacing, "lines"), strip.placement = "outside",  
      strip.text.y = element_text(angle = 0),
      strip.text.x = element_text(colour = "white", face = "bold"), 
      strip.background = element_rect(colour = "white")
    )
  plot_strips <- plot_adjusted_strips(plot, strip_colors = c(
    "#2EBAED", "#000000", "#DE1C14",
    "#D4D2D2", "#ADCC54", "#F0D0CE"
  ))
  return(plot_strips)
}



# ggsave doesnt export the facet strip color properly
plot_indel_contexts2 <- function (counts, same_y = FALSE, extra_labels = FALSE, condensed = FALSE, striplabel = "top", stripcolor = "fill",
                                  basesize = 6, legend = FALSE, show_indel_counts_strips = TRUE, y_label = " (Fraction)") 
{
  count <- muttype <- muttype_sub <- muttype_total <- sample <- NULL
  counts <- counts %>% as.data.frame() %>% tibble::rownames_to_column("muttype_total") %>% 
    tidyr::separate(muttype_total, c("muttype", "muttype_sub"), 
                    sep = "_(?=[0-9])") %>% dplyr::mutate(muttype = factor(muttype, 
                                                                           levels = unique(muttype))) %>% tidyr::gather(key = "sample", 
                                                                                                                        value = "count", -muttype, -muttype_sub) %>% dplyr::mutate(sample = factor(sample, levels = unique(sample)))
  counts$Label <- ifelse(grepl("deletion", x = counts$muttype) == T, "del", "ins")
  
  # Add the categories to the indel counts (these will used as the facet strip headers)
  indel_categories <- data.frame(muttype = levels(counts$muttype),
                                 category = c(rep("1bp Deletion", 2),
                                              rep("1bp Insertion", 2),
                                              rep(">1bp Deletion at repeats\n(Deletion length)", 4),
                                              rep(">1bp Insertion at repeats\n(Insertion length)", 4),
                                              rep("Deletion with MH\n(Deletion length", 4)),
                                 xaxis = c(rep("Homopolymer length", 2),
                                           rep("Homopolymer length", 2),
                                           rep("Number of repeat units", 4),
                                           rep("Number of repeat units", 4),
                                           rep("Microhomology\nlength", 4)))
  indel_categories$category <- factor(indel_categories$category, levels = c("1bp Deletion", "1bp Insertion",">1bp Deletion at repeats\n(Deletion length)",
                                                                            ">1bp Insertion at repeats\n(Insertion length)",
                                                                            "Deletion with MH\n(Deletion length"))
  indel_categories$xaxis <- factor(indel_categories$xaxis, levels = c("Homopolymer length", "Number of repeat units",  "Microhomology\nlength"))
  
  counts <- merge(counts,indel_categories, by = "muttype", all.x = T)
  
  nr_muts <- counts %>% dplyr::group_by(sample) %>% dplyr::summarise(nr_muts = round(sum(count)))
  if(show_indel_counts_strips == TRUE){
    facet_labs_y <- stringr::str_c(nr_muts$sample, "\n(n = ", 
                                   nr_muts$nr_muts, ")")
  } else {
    facet_labs_y <- stringr::str_c(nr_muts$sample)
  }
  
  names(facet_labs_y) <- nr_muts$sample
  facet_labs_x <- c("C", "T", "C", "T", 2, 3, 4, 
                    "5+", 2, 3, 4, "5+", 2, 3, 4, "5+")
  names(facet_labs_x) <- levels(counts$muttype)
  if (same_y) {
    facet_scale <- "free_x"
  }
  else {
    facet_scale <- "free"
  }
  
  
  indel_colors <- c( "#FDBE6F", "#FF8001", "#B0DD8B", "#36A12E", "#FDCAB5", 
                     "#FC8A6A", "#F14432", "#BC141A", "#D0E1F2", "#94C4DF", 
                     "#4A98C9", "#1764AB", "#E2E2EF", "#B6B6D8", "#8683BD", 
                     "#61409B")
  strip_colors <-  c(rep("gray90", 8), "#FDBE6F", "#FF8001", "#B0DD8B", "#36A12E", "#FDCAB5", 
                     "#FC8A6A", "#F14432", "#BC141A", "#D0E1F2", "#94C4DF", 
                     "#4A98C9", "#1764AB", "#E2E2EF", "#B6B6D8", "#8683BD", 
                     "#61409B")
  if (extra_labels) {
    title <- stringr::str_c("Deletion           ", "Insertion          ", 
                            "Deletion                                   ", "Insertion                                  ", 
                            "Deletion (MH)")
    x_lab <- stringr::str_c("Homopolymer length                            ", 
                            "Number of repeat units                                                                               ", 
                            "Microhomology length")
  }
  else {
    title <- x_lab <- ""
  }
  if (condensed == TRUE) {
    width <- 1
    spacing <- 0
  }
  else {
    width <- 0.6
    spacing <- 0.5
  }
  
  if(legend == FALSE){
    legend_position <- "none"
  } else {
    legend_position <- "right"
  }
  
  if(striplabel == "top"){
    fig <- ggplot(counts, aes(x = muttype_sub, y = count, fill = muttype, 
                              width = width)) + 
      geom_bar(stat = "identity", col = "white", size = 0.2) + 
      facet_nested(sample ~ xaxis +category + muttype, scales = facet_scale, space = "free_x", labeller = labeller(muttype = facet_labs_x, sample = facet_labs_y) ) + 
      scale_fill_manual(values = indel_colors) + 
      scale_y_continuous(expand = expansion(mult = c(0,0.1))) +   
      theme_classic(base_size = basesize) + 
      labs(fill = "Mutation type", title = title, y = paste("Indels", y_label, sep = ""), x = x_lab) +
      annotation_custom(grid::linesGrob(x = c(0, 0), gp = grid::gpar(lwd = 0.5, col = "lightgrey"))) +
      annotation_custom(grid::linesGrob(y = c(0, 0), gp = grid::gpar(lwd = 1.5, col = "black"))) +
      theme(panel.grid.major.x = element_blank(),
            panel.grid.minor.y = element_blank(), 
            panel.grid = element_blank(),
            panel.spacing.x = unit(spacing, "lines"),
            #panel.border = element_rect(color = "lightgrey", fill = NA),
            strip.background.x = element_rect(color = "white"),
            strip.background.y = element_rect(fill = "gray92", color = "white"),
            axis.text.x = element_text(size = 4),
            strip.placement = "outside", 
            strip.text.y.right = element_text(angle = 0),
            legend.position = legend_position) 
    
  } else {
    fig <- ggplot(counts, aes(x = muttype_sub, y = count, fill = muttype, 
                              width = width)) + geom_bar(stat = "identity") + 
      facet_nested(sample ~ category +muttype, scales = facet_scale, space = "free_x", labeller = labeller(muttype = facet_labs_x, sample = facet_labs_y), switch = "x") + 
      cale_fill_manual(values = indel_colors) + 
      scale_y_continuous(expand = expansion(mult = c(0,0.1))) +   
      theme_bw(base_size = basesize) + 
      labs(fill = "Mutation type", title = title, y = paste("Indels", y_label, sep = ""), x = x_lab) +
      theme(panel.grid.major.x = element_blank(),
            panel.grid.minor.y = element_blank(), 
            panel.grid = element_blank(),
            panel.spacing.x = unit(spacing, "lines"), 
            strip.background = element_rect(color = "white"),
            strip.placement = "outside", 
            strip.text.y.right = element_text(angle = 0),
            legend.position = legend_position)
  }
  g <- ggplot_gtable(ggplot_build(fig))
  if(striplabel == "top"){
    strip_both <- which(grepl('strip-t', g$layout$name))
  } else {
    strip_both <- which(grepl('strip-b', g$layout$name))
  }
  if(stripcolor == "fill"){
    fills <- strip_colors
  } else{
    fills <- ifelse(grepl("deletion", x = levels(counts$muttype)) == T, "#FB9A99", "#A6CEE3")
  }
  k <- 1
  
  for (i in strip_both) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  plot(g)
  ## Return doesn't work with the edited facet strips
  #return(fig)
}



# This function reads a vcf and selects indels, including the VAF
readIndels <- function(vcf, chromosomes = c(1:22, "X", "Y"), rename_indels = F, genome = "hg38", VAF_threshold = 0.2){
  variants <- readVcf(vcf, genome = genome)
  indels <- variants[isIndel(variants),]
  indels <- indels[which(as.vector(seqnames(indels)) %in% chromosomes),]
  indels_gr <- VCF_to_GR(indels, chromosomes = chromosomes)
  indels_gr$VAF <- as.numeric(geno(indels)$VAF)
  
  # Some variants may have SNP ids as name, instead of chr-pos-indel
  if(rename_indels == TRUE){
    names(indels_gr) <- paste(gsub(pattern = "chr", replacement = "", as.vector(seqnames(indels_gr))), ":", start(indels_gr), "_",  as.vector(indels_gr$REF), "/",as.vector(unlist(indels_gr$ALT)), sep = "")
  }
  indels_gr_filtered <- indels_gr[which(indels_gr$VAF >= VAF_threshold),]
  print(paste("# Removed ", length(indels_gr)-length(indels_gr_filtered), " indels with VAF <", VAF_threshold, sep = ""))
  return(indels_gr_filtered)
}



### Function to filter indels present in the blacklist and in homopolymers. 
filter_indels <- function(indel_grl, indel_recurrency_grl, remove_homopolymer = TRUE, mode = "filter"){
  indels_filtered <- list()
  for(sample in names(indel_grl)){
    print(sample)
    indels_sample <- indel_grl[[sample]]
    indels_sample$repeat_length <-  (width(unlist(indels_sample$ALT))-1) * indels_sample$muttype_sub
    print(paste("### Input ", length(indels_sample), " indels", sep = ""))
    
    olap_blacklist <- findOverlaps(indels_sample, indel_recurrency_grl)
    print(paste("### Removed ", length(olap_blacklist), " indels overlapping with Blacklist", sep = ""))
    if(length(olap_blacklist) > 0 & mode == "filter"){
      indels_filtered_sample <- indels_sample[-queryHits(olap_blacklist)]
    } else if (length(olap_blacklist) > 0 & mode == "flag") {
      indels_filtered_sample <- indels_sample
      indels_filtered_sample$FILTER[queryHits(olap_blacklist)] <- "RECURRENT"
    } else {
      indels_filtered_sample <- indels_sample
    }
    ### The current blacklist does not remove all 1bp and 2bp insertions in homopolymers/repeats. Remove these as well
    if(remove_homopolymer == TRUE){
      remove <- which(indels_filtered_sample$repeat_length > 4 & grepl(indels_filtered_sample$muttype, pattern = "insertion"))
      print(paste("### Removed ", length(remove), " insertions in homopolymers"))
      if(length(remove) > 0 & mode == "filter"){
        indels_filtered_sample <- indels_filtered_sample[-remove,]
      } else if (length(remove) > 0 & mode == "flag") {
        indels_filtered_sample$FILTER[remove][indels_filtered_sample$FILTER[remove] != "PASS"] <- paste(indels_filtered_sample$FILTER[remove][indels_filtered_sample$FILTER[remove] != "PASS"], "HOMOPOLYMER", sep = ";")
        indels_filtered_sample$FILTER[remove][indels_filtered_sample$FILTER[remove] == "PASS"] <- "HOMOPOLYMER"
      }
    }
    print(paste("### Output ", length(indels_filtered_sample), " indels", sep = ""))
    indels_filtered[[sample]] <- indels_filtered_sample
  }
  return(indels_filtered)
}



# Standard set of colors used 
walker_colors <- c("#D55E00", "#669bbc", "#009e73", "gray")
names(walker_colors) <- c("Low", "Intermediate", "High", NA)

# This function reads the Walker output vcf.gz files in the "walker_vcf_directory". "samples" should be in the vcf.gz filenames.
# Can only run for snv or indel seperately
# Output is a dataframe with all the variants with their walker scores
# Thresholds are now determine arbitarily
RetrieveWalkerScores <- function(walker_vcf_directory, samples, type = "snv",
                                 high_score_threshold = 1000, 
                                 low_score_threshold = 1){
  walker_vcf_files <- list.files(walker_vcf_directory, pattern = ".vcf.gz$", full.names = T)
  walker_scores <- data.frame()
  for(sample in samples){
    print(sample)
    walker_vcf_file <- walker_vcf_files[grep(sample, walker_vcf_files)]
    print(walker_vcf_file)
    # Read the Walker VCF
    walker_vcf <- readVcf(walker_vcf_file, genome = "hg38")
    # Select only the snv or the indel
    if(type == "snv"){
      walker_variants <- walker_vcf[isSNV(walker_vcf),]
    } else if (type == "indel"){
      walker_variants <- walker_vcf[isIndel(walker_vcf),]
    }
    # Select autosomes only
    walker_variants <- walker_variants[as.vector(seqnames(walker_variants)) %in% c(1:22),]
    
    # Collect all the walker scores
    walker_scores_sample <- data.frame(Variant = row.names(geno(walker_variants)$WS), 
                                       WS = as.numeric(geno(walker_variants)$WS), check.names = F)
    walker_scores_sample$WS_Count <- walker_scores_sample$WS
    walker_scores_sample$WS_Count[which(as.numeric(walker_scores_sample$WS) < low_score_threshold)] <- "Low"
    walker_scores_sample$WS_Count[which(as.numeric(walker_scores_sample$WS) > high_score_threshold)] <- "High"
    walker_scores_sample$WS_Count[which(as.numeric(walker_scores_sample$WS) >= low_score_threshold  &as.numeric(walker_scores_sample$WS) <= high_score_threshold)] <- "Intermediate"
    walker_scores_sample$WS_Count[is.na(walker_scores_sample$WS_Count)] <- NA
    walker_scores_sample$WS_Count <- factor(walker_scores_sample$WS_Count, levels = c(NA, "Low", "Intermediate", "High"))
    walker_scores_sample$sample <- sample
    walker_scores <- rbind(walker_scores, walker_scores_sample)
  }
  return(walker_scores)
}

# This calculaes the number of variants for each category (Low, Intermediate, High, NA) per sample
CalculateWalkerCounts <- function(walker_scores){
  walker_counts <- walker_scores %>% group_by(sample, WS_Count) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n))
  return(walker_counts)
}


# The Walker scores are extrapolated by multiplying the amount of undetermined variants (WS_Count=NA) by the frequency of the called variants (WS_Count!=NA)
Extrapolate_Walker <- function(walker_counts, walker_scores){
  walker_counts_extrapolated <- data.frame()
  walker_counts_called <- walker_scores[!is.na(walker_scores$WS_Count),] %>% group_by(sample, WS_Count) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n))
  
  ## Extrapolate walker
  for(sample in unique(walker_counts$sample)){
    #print(sample)
    walker_counts_uncalled_sample <- walker_counts[walker_counts$sample == sample & is.na(walker_counts$WS_Count),]
    walker_counts_called_sample <- walker_counts_called[walker_counts_called$sample == sample,]
    walker_counts_extrapolated_sample <- data.frame(sample = sample, 
                                                    WS_Count = walker_counts_called_sample$WS_Count, 
                                                    n = walker_counts_uncalled_sample$n * walker_counts_called_sample$freq)
    walker_counts_extrapolated_sample$Extrapolated <- walker_counts_extrapolated_sample$n + walker_counts_called_sample$n
    #print(walker_counts_extrapolated_sample)
    walker_counts_extrapolated <- rbind(walker_counts_extrapolated, walker_counts_extrapolated_sample)
  }
  return(walker_counts_extrapolated)
}


plot_walker_counts <- function(walker_counts, mode = "absolute"){
  if(mode == "absolute"){
    plot <- ggplot(walker_counts, aes(x = sample, fill = WS_Count, y = n)) + 
      geom_bar(stat = "identity") + 
      theme_classic(base_size = 8) +
      theme(axis.text.x = element_text(angle = 60,  hjust=1)) + 
      scale_y_continuous(expand = c(0, 0)) + 
      scale_fill_manual(values = walker_colors) + 
      labs(y = "Somatic variants (Autosomal)", fill = "Walker Score", x = "Sample")
  } else {
    plot <- ggplot(walker_counts, aes(x = sample, fill = WS_Count, y = freq)) + geom_bar(stat = "identity") +
      geom_text(aes(label = paste(round(freq*100,0), "%", sep =""), y=freq), position = position_stack(vjust = 0.5), size= 3) +
      theme_classic(base_size = 8) + 
      scale_y_continuous(expand = c(0, 0)) + 
      scale_fill_manual(values = walker_colors) +
      theme(axis.text.x = element_text(angle = 60,  hjust=1)) + 
      labs(y = "Somatic variants (% of total)", fill = "Walker Score", x = "Sample")
  }
  return(plot)
}

plot_Walker_Extrapolated <- function(Walker_Extrapolated){
  plot <- ggplot(Walker_Extrapolated, aes(x = sample, fill = WS_Count, y = Extrapolated)) + 
    geom_bar(stat = "identity", col = "black") +
    geom_text(aes(label = round(Extrapolated,0), y = Extrapolated), position = position_stack(vjust = 0.5), size= 3) +
    theme_classic(base_size = 8) + 
    scale_y_continuous(expand = c(0, 0)) + 
    scale_fill_manual(values = walker_colors) +
    theme(axis.text.x = element_text(angle = 60,  hjust=1)) + 
    labs(y = "Somatic variants (Extrapolated autosomes)", fill = "Walker Score", x = "Sample")
  return(plot)
}

plot_Walker_distribution <- function(walker_scores){
  plot <- ggplot(walker_scores, aes(x = sample, y = WS, fill = sample)) +
    geom_hline(yintercept = c(100, 1000), linetype = 2, col = "grey") +
    geom_boxplot() +  
    coord_cartesian(ylim = c(0,3000)) + 
    theme_classic(base_size = 8) + 
    theme(axis.text.x = element_text(angle = 60,  hjust=1), legend.position = "none") +
    labs(y = "Walker Score", x = "Sample")
  return(plot)
}


# Makes a list with two mutation matrices. Each matrix contains the mutation profiles at different PTAprobsCutoffs. One matrix contains the mutations passing the filter (eg PTAprob of variant <= Cutoff), the other contains the mutations failing the filter (eg PTAprob of variant > Cutoff)
# Additionally, the original PTAprobCutoff as determined by PTATO is also included in the list
calc_cosim_mutmat <- function(gr, cutoffs =  seq(0.1, 0.8, 0.025)){
  grl_cosim <- list()
  grl_removed <- list()
  
  for(cutoff in cutoffs){
    #print(cutoff)
    grl_cosim[[as.character(cutoff)]] <- gr[which(gr$PTAprob <= cutoff)]
    grl_removed[[as.character(cutoff)]] <- gr[which(gr$PTAprob > cutoff)]
  } 
  mut_mat <- MutationalPatterns::mut_matrix(vcf_list = grl_cosim, ref_genome = ref_genome)
  mut_mat_removed <- MutationalPatterns::mut_matrix(vcf_list = grl_removed, ref_genome = ref_genome)
  mut_mat_sample <- list(PASS = mut_mat, FAIL = mut_mat_removed, PTAprobsCutoff = unique(gr$PTAprobCutoff))
  return(mut_mat_sample)
}
