### Figure S7
library(VariantAnnotation)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(stringr)
library(MutationalPatterns)
library(dplyr)

ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

Input_dir <- "/path/to/MendeleyData/directory/"
Scripts_dir <-  "/path/to/GithubScripts/directory/"
Output_dir<-  "/path/to/output/directory/"
Resources_dir <- "/path/to/MendeleyData/directory/Resources/
"
if(dir.exists(Output_dir) == F) { dir.create(Output_dir)}
source(paste(Scripts_dir, "GeneralFunctions.R", sep = ""))


### Functions
get_sample_names = function(vcf_fnames, 
                            str_remove = ".ptato.*"){
  sample_names = vcf_fnames %>% 
    basename() %>% 
    str_remove(str_remove)
  sample_names <- gsub(".*_", "", sample_names)
  return(sample_names)
}

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

# Determine file names
snv_unfilt_vcfs_fnames = list.files(file.path(paste(Input_dir, "GonzalezPena2021_CordBlood/PTATO/", sep = ""), "snvs", "Cord"), pattern = ".vcf.gz$", full.names = TRUE)
snv_unfilt_sample_names = get_sample_names(snv_unfilt_vcfs_fnames, str_remove = ".snvs.ptato.*")
names(snv_unfilt_vcfs_fnames) <- snv_unfilt_sample_names
snv_filt_vcfs_fnames = list.files(file.path(paste(Input_dir, "GonzalezPena2021_CordBlood/PTATO/", sep = ""), "snvs", "Cord"), pattern = "filtered.vcf.gz$", full.names = TRUE, recursive = TRUE)
snv_filt_sample_names = get_sample_names(snv_filt_vcfs_fnames, str_remove = ".snvs.ptato.*")
names(snv_filt_vcfs_fnames) <- snv_filt_sample_names

grl_raw <- list()
for(sample in snv_filt_sample_names){
  print(sample)
  
  grl_raw[[sample]] <- readPTATOvcf(vcf_unfiltered = snv_unfilt_vcfs_fnames[sample],
                              vcf_filtered = snv_filt_vcfs_fnames[sample], 
                              VAF_threshold = 0.2,
                              type = "snv", chromosomes = c(1:22))
}

grl_raw <- lapply(grl_raw, function(x) x[which(x$FILTER != "FAIL_VAF"),])
grl_PASS <- lapply(grl_raw, function(x) x[which(x$FILTER == "PASS"),])
grl_FAIL <- lapply(grl_raw, function(x) x[which(x$FILTER == "FAIL"),])

### Read the SCAN2 output
SCAN2_rescuedmuts <- read.delim(paste(Input_dir, "GonzalezPena2021_CordBlood/SCAN2/panel/run_job/rescue-alll/rescue/rescued_muts.txt", sep = ""))

# This function changes a SCAN2 rescue_muts.txt dataframe to granges for MutationalPatterns
SCAN2_to_GRanges <- function(rescued_muts, type = ""){
  SCAN2_grl <- list()
  for(sample in unique(rescued_muts$sample)){
    print(sample)
    rescued_muts_sample <- rescued_muts[which(rescued_muts$sample == sample),]
    SCAN2_g <- GRanges(seqnames = rescued_muts_sample$chr, 
                       IRanges(start = rescued_muts_sample$pos, 
                               end = rescued_muts_sample$pos, 
                               strand = "*"), 
                       REF = rescued_muts_sample$refnt, 
                       ALT = rescued_muts_sample$altnt,
                       QUAL = NA,
                       FILTER = ifelse(rescued_muts_sample$pass == TRUE, "PASS", "RESCUE"),
                       VAF = rescued_muts_sample$af)
    names(SCAN2_g) <- paste(gsub(pattern = "chr", replacement = "", rescued_muts_sample$chr), ":", 
                            rescued_muts_sample$pos, "_", 
                            rescued_muts_sample$refnt, "/",
                            rescued_muts_sample$altnt, sep = "")
    genome(SCAN2_g) <- "hg38"
    if(type == "snv"){
      SCAN2_g <- get_mut_type(SCAN2_g, type = "snv")
    }
    SCAN2_grl[[sample]] <- SCAN2_g
  }
  return(SCAN2_grl)
}

sSNVs_SCAN2 <- SCAN2_to_GRanges(rescued_muts = SCAN2_rescuedmuts, type = "snv")


snv_counts <- data.frame(Sample = names(grl_raw), Unfiltered = lengths(grl_raw), PTATO = lengths(grl_PASS), Removed = lengths(grl_FAIL))

snv_counts$Group <- str_remove(snv_counts$Sample, "-[0-9]$")
snv_counts$Treatment <- str_remove(snv_counts$Group, "-[0-9]$")
snv_counts$Dose <- str_remove(snv_counts$Group, ".*-")
snv_counts$Dose_label <- ifelse(snv_counts$Dose == 0, "None", "Low" )
snv_counts$Dose_label[snv_counts$Dose == 3 & snv_counts$Treatment == "MAN"] <- "Moderate"
snv_counts$Dose_label[snv_counts$Dose == 2 & snv_counts$Treatment == "ENU"] <- "Moderate"

snv_counts$Dose_label[snv_counts$Dose == 3 & snv_counts$Treatment == "ENU"] <- "High"

snv_counts_scan2 <- data.frame(Sample = names(sSNVs_SCAN2), SCAN2 = lengths(sSNVs_SCAN2))
snv_counts_scan2$Group <- str_remove(snv_counts_scan2$Sample, "-[0-9]$")
snv_counts_scan2$Treatment <- str_remove(snv_counts_scan2$Group, "-[0-9]$")

snv_counts_merged <- merge(snv_counts, snv_counts_scan2, by = c("Sample", "Group", "Treatment"), all.x = T)

snv_counts_m <- melt(snv_counts_merged)
snv_counts_m$Treatment <- factor(snv_counts_m$Treatment, levels = c("VHC", "MAN", "ENU"))
snv_counts_m$Dose_label <- factor(snv_counts_m$Dose_label, levels = c("None", "Low", "Moderate", "High"))
snv_counts_m$variable <- factor(snv_counts_m$variable, levels = c("Unfiltered", "Removed", "PTATO", "SCAN2"))


snv_counts_means <- snv_counts_m %>% dplyr::group_by(Treatment, Dose, variable, Dose_label) %>% dplyr::summarise(mean = round(mean(value),0), 
                                                                                                                 y_pos = max(value)+300)

ClassificationColors <- c( "#377eb8","#e41a1c", "#4daf4a", "orange")
names(ClassificationColors) <- c( "Unfiltered", "Removed","PTATO", "SCAN2")

Treatment_colors <- c("#66c2a5", "#fc8d62", "#8da0cb")
names(Treatment_colors) <- c("VHC", "MAN", "ENU")

SNV_counts_plots <- ggplot(snv_counts_m, aes(x = Dose_label, y = value, col = variable, fill = variable)) +
  geom_boxplot(outlier.shape =  NA, size = 0.3, position = position_dodge(preserve = "single")) +
  geom_point(aes( fill = variable),pch = 21, position = position_jitterdodge(jitter.width = 0.2), col = "black", size = 0.8) +
  geom_text(data = snv_counts_means, aes(x = Dose_label, y = y_pos, label = mean, col = variable), 
            position = position_dodge(width = 0.8), size = 1.6, fontface = "italic") +
  facet_grid(.~Treatment, space = "free", scale = "free", switch = "x") +
  scale_color_manual(values = ClassificationColors) +
  scale_fill_manual(values = alpha(ClassificationColors, 0.1)) + 
  theme_classic(base_size = 6) +
  theme(strip.placement = "outside", strip.background = element_rect( color = "white"), strip.text = element_text(color = "white", face = "bold")) +
  labs(x = "Treatment dose", y = "Single Base Substitutions\nper sample", col = "Filtering", fill = "Filtering")

pdf(file = paste(Output_dir, "FigureS7A_SNV_Counts.pdf",sep = ""), width = 4.8, height = 2, onefile = F)

plot_adjusted_strips(SNV_counts_plots, strip_colors = c(
  Treatment_colors
))
dev.off()


## Signatures
signatures = get_known_signatures()
cosmic_sig_fname = paste(Resources_dir, "sigProfiler_SBS_working_signatures_incl_hspc.txt", sep = "")
cosmic_sign = read.table(cosmic_sig_fname, sep = "\t", header = T)
cosmic_sign = as.matrix(cosmic_sign[,-c(1,2), drop = FALSE])
HSPC <- cosmic_sign[,"HSPC", drop = FALSE]
pta_sig_fname = paste(Resources_dir, "PTA_Artefact_Signature.txt", sep = "")
pta_sign = read.table(pta_sig_fname, sep = "\t", header = T)
PTA = as.matrix(pta_sign[,-1, drop = FALSE])
sig_enu = get_known_signatures(source = "SIGNAL", sig_type = "exposure")[,18,drop=F]
colnames(sig_enu) = "ENU"
signatures <- do.call(cbind, list(signatures, PTA, HSPC, sig_enu))

### 96-profiles
unfilt_mut_mat = mut_matrix(grl_raw, BSgenome.Hsapiens.UCSC.hg38)
unfilt_mut_mat_group = pool_mut_mat(unfilt_mut_mat, snv_counts$Group )
plot_96_profile(unfilt_mut_mat_group, condensed = TRUE, ymax = 0.1)


filt_mut_mat = mut_matrix(grl_PASS, BSgenome.Hsapiens.UCSC.hg38)
filt_mut_mat_group = pool_mut_mat(filt_mut_mat, snv_counts$Group )
plot_96_profile(filt_mut_mat_group, condensed = TRUE, ymax = 0.1)

removed_mut_mat <- mut_matrix(grl_FAIL, BSgenome.Hsapiens.UCSC.hg38)
removed_mut_mat_group = pool_mut_mat(removed_mut_mat, snv_counts$Group )

scan2_mutmat <-  mut_matrix(sSNVs_SCAN2, BSgenome.Hsapiens.UCSC.hg38)
scan2_mutmat_group = pool_mut_mat(scan2_mutmat, snv_counts$Group[snv_counts$Group %in% c("ENU-3", "VHC-0")] )


mut_mat_comp = cbind(unfilt_mut_mat_group[,c("VHC-0"), drop = FALSE], 
                     filt_mut_mat_group[,c("VHC-0"), drop = FALSE], 
                     removed_mut_mat_group[,c("VHC-0"), drop = FALSE],
                     scan2_mutmat_group[,c("VHC-0"), drop = FALSE],
                     
                     unfilt_mut_mat_group[,c("MAN-3"), drop = FALSE],
                     filt_mut_mat_group[,c("MAN-3"), drop = FALSE],
                     removed_mut_mat_group[,c("MAN-3"), drop = FALSE],
                     
                     unfilt_mut_mat_group[,c("ENU-3"), drop = FALSE],
                     filt_mut_mat_group[,c("ENU-3"), drop = FALSE],
                     removed_mut_mat_group[,c("ENU-3"), drop = FALSE],
                     scan2_mutmat_group[,c("ENU-3"), drop = FALSE])
colnames(mut_mat_comp) <- c("VHC-none\nUnfiltered", "VHC-none\nPTATO","VHC-none\nRemoved", "VHC-none\nSCAN2",
                            "MAN-moderate\nUnfiltered", "MAN-moderate\nPTATO","MAN-moderate\nRemoved",
                            "ENU-high\nUnfiltered", "ENU-high\nPTATO", "ENU-high\nRemoved", "ENU-high\nSCAN2")

mut_mat_comp_withENU <- cbind(mut_mat_comp, signatures[,"ENU"])
colnames(mut_mat_comp_withENU)[ncol(mut_mat_comp_withENU)] <- "ENU"


### 96-profiles
strip_colors_row <- ifelse(grepl("ENU", colnames(mut_mat_comp_withENU)), Treatment_colors["ENU"], Treatment_colors["VHC"])
strip_colors_row[grep("MAN", colnames(mut_mat_comp_withENU))] <- Treatment_colors["MAN"]

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
    "#D4D2D2", "#ADCC54", "#F0D0CE", strip_colors_row
  ))
  return(plot_strips)
}

pdf(file = paste(Output_dir, "FigureS7B_96profiles.pdf",sep = ""), width = 3.5, height = 4.5, onefile = F)
plot_96_profile_small(mut_mat_comp_withENU, ymax = 0.1, condensed = T)
dev.off()

### Cosine similarities
cos_sim_samples_signatures <- cos_sim_matrix(mut_mat_comp, signatures[,c("PTA", "SBS5", "ENU")])

cos_sim_samples_signatures_m <- melt(cos_sim_samples_signatures)
cos_sim_samples_signatures_m$Filter <- gsub(cos_sim_samples_signatures_m$Var1, pattern = ".*\n", replacement = "")
cos_sim_samples_signatures_m$Sample <- gsub(cos_sim_samples_signatures_m$Var1, pattern = "\n.*", replacement = "")

cos_sim_samples_signatures_m$Sample <-  factor(cos_sim_samples_signatures_m$Sample, levels = c("VHC-none",  "MAN-moderate", "ENU-high" ))
cos_sim_samples_signatures_m$Var2 <- factor(cos_sim_samples_signatures_m$Var2, levels = c("SBS5", "ENU", "PTA"))
cos_sim_samples_signatures_m$Filter <- factor(cos_sim_samples_signatures_m$Filter, levels = c("SCAN2", "PTATO","Removed" ,"Unfiltered"))

cossim_plot <- ggplot(cos_sim_samples_signatures_m, aes(x = Var2, y = Filter, fill = value, label = round(value,2))) + geom_raster() +
  geom_text(size = 1.8) + 
  facet_grid(Sample~., switch = "both", scales = "free") +
  scale_fill_gradient(low = "white", high = "firebrick3") +
  scale_x_discrete(position = "top", expand = c(0,0)) +
  theme_classic(base_size = 6) +
  theme( strip.placement = "outside", axis.text.x = element_text(angle = -45, hjust = 1), 
         strip.background = element_rect(color = "white", fill = "lightgray")) +
  labs(x = "Mutational Profile", y= "Sample", fill = "Cosine\nSimilarity") + 
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 5))

pdf(file = paste(Output_dir, "FigureS7C_CosineSimilarity_Heatmap.pdf", sep = ""), width = 2.2, height = 2.5, onefile = F)
plot_adjusted_strips(cossim_plot, strip_colors = c(Treatment_colors))
dev.off()

### Bootstrapped refitting

makeBootstrapMatrix <- function(mutmat, signatures, n_boots = 100, add_names = T){
  contri_boots <- fit_to_signatures_bootstrapped(mutmat,
                                signatures,
                                n_boots = 100,
                                method = "strict", max_delta = 0.005)
  # Change to long format for plotting
  contri_tb <- contri_boots %>% as.data.frame() %>% tibble::rownames_to_column("exp") %>% 
    tidyr::gather(key = "sig", value = "contri", -exp) %>% 
    dplyr::mutate(Sample = gsub("_[^_]+$", "", exp), Sample = factor(Sample, 
                                                                     levels = unique(Sample)), sig = factor(sig, levels = unique(sig)))
  # Calculate the mean contribution for each signature
  contri_tb2 <- contri_tb %>% dplyr::group_by(Sample, sig) %>% 
    dplyr::summarise(mean = mean(contri)) %>% 
    dplyr::ungroup()
  
  if(add_names == TRUE){
    contri_tb2$Treatment <- str_remove(contri_tb2$Sample, "-[0-9]$")
    contri_tb2$Dose <- str_remove(contri_tb2$Sample, ".*-")
    contri_tb2$Dose_label <- ifelse(contri_tb2$Dose == 0, "None", "Low" )
    contri_tb2$Dose_label[contri_tb2$Dose == 3 & contri_tb2$Treatment == "MAN"] <- "Moderate"
    contri_tb2$Dose_label[contri_tb2$Dose == 2 & contri_tb2$Treatment == "ENU"] <- "Moderate"
    
    contri_tb2$Dose_label[contri_tb2$Dose == 3 & contri_tb2$Treatment == "ENU"] <- "High"
    
  }
  
  return(contri_tb2)
}

contri_boots_unfiltered <- makeBootstrapMatrix(mutmat = unfilt_mut_mat_group, signatures = signatures[,c("PTA", "SBS1", "SBS5", "ENU")])
contri_boots_unfiltered$Label <- "Unfiltered"

contri_boots_filtered <- makeBootstrapMatrix(mutmat = filt_mut_mat_group, signatures = signatures[,c("PTA", "SBS1", "SBS5", "ENU")])
contri_boots_filtered$Label <- "PTATO"

contri_boots_removed <- makeBootstrapMatrix(mutmat = removed_mut_mat_group, signatures = signatures[,c("PTA", "SBS1", "SBS5", "ENU")])
contri_boots_removed$Label <- "Removed"

contri_boots_scan2 <- makeBootstrapMatrix(mutmat = scan2_mutmat_group, signatures = signatures[,c("PTA", "SBS1", "SBS5", "ENU")])
contri_boots_scan2$Label <- "SCAN2"

contri_boots_merged <- rbind(contri_boots_unfiltered, contri_boots_removed, contri_boots_filtered, contri_boots_scan2)
contri_boots_merged$Facet <- paste(contri_boots_merged$Treatment, contri_boots_merged$Dose_label, sep = "\n")


contri_boots_merged$Facet <- factor(contri_boots_merged$Facet, levels = c("VHC\nNone",  "MAN\nLow","MAN\nModerate", "ENU\nLow", "ENU\nModerate","ENU\nHigh"))
contri_boots_merged$Label <- factor(contri_boots_merged$Label, levels = c("Unfiltered", "Removed","PTATO", "SCAN2"))

signature_colors <- c("#E7298A","#7570B3","#66A61E","#D95F02" )
names(signature_colors) <- c("PTA", "SBS1", "SBS5", "ENU")


refit_plot <- ggplot(contri_boots_merged, aes( x = Label, y= mean)) + 
  geom_bar(stat = "identity", aes(fill = sig), col = "black", linewidth = 0.2) +
  facet_grid(.~Facet, switch = "x", scales = "free", space = "free") + 
  # geom_text(aes(label = Label, fill = Signature), position = position_stack(vjust = 0.5), size = 1.8, col = "white") +
  # geom_text(data = strict_res_sum, aes(y = Label + 200, label = Label), size = 1.8) +
  theme_classic(base_size = 6) +
  scale_y_continuous(expand = c(0,0), limits = c(0,22000)) + 
  scale_fill_manual(values = signature_colors) +
  theme(strip.placement = "outside", 
        strip.background = element_rect(color = "white", fill = "lightgray"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.key.size = unit(0.25, 'cm')) +
  labs(x = "Treatment dose", y = "Absolute Contribution\n(Strict Refit)", fill = "Signature")


pdf(file = paste(Output_dir, "FigureS7D_StrictRefit.pdf", sep = ""), width = 3.6, height = 2, onefile = F)
plot_adjusted_strips(refit_plot, strip_colors = c(Treatment_colors[c(1,2,2,3,3,3)]
))
dev.off()

contri_boots_filtered$Treatment_label <- paste(contri_boots_filtered$Treatment, contri_boots_filtered$Dose_label, sep = "\n")

PrecRecall <- data.frame()  
# Sensitivity
for(Sample in unique(contri_boots_filtered$Sample)){
  print(Sample)
  
  # Precision (how many variants after filtering are PTA artefacts?)
  # Sensitivity (how many non-PTA variants detected in the unfiltered samples are also found in the filtered samples?)
  SampleDF <- data.frame(Sample = Sample, 
                         Treatment = contri_boots_filtered$Treatment_label[contri_boots_filtered$Sample == Sample],
                         Precision = 1 - sum(contri_boots_filtered$mean[contri_boots_filtered$sig == "PTA" & contri_boots_filtered$Sample == Sample]) /
                           sum(contri_boots_filtered$mean[contri_boots_filtered$Sample == Sample]),
                         
                         Sensitivity = sum(contri_boots_filtered$mean[contri_boots_filtered$sig != "PTA" & contri_boots_filtered$Sample == Sample]) / 
                           sum(contri_boots_unfiltered$mean[contri_boots_unfiltered$sig != "PTA" & contri_boots_unfiltered$Sample == Sample]),
                         Caller = "PTATO")
  
  PrecRecall <- rbind(PrecRecall, SampleDF)
  
  SampleDF_SCAN2 <- data.frame(Sample = Sample, 
                               Treatment = contri_boots_filtered$Treatment_label[contri_boots_filtered$Sample == Sample],
                
                         Precision = 1 - sum(contri_boots_scan2$mean[contri_boots_scan2$sig == "PTA" & contri_boots_scan2$Sample == Sample]) /
                           sum(contri_boots_scan2$mean[contri_boots_scan2$Sample == Sample]),
                         
                         Sensitivity = sum(contri_boots_scan2$mean[contri_boots_scan2$sig != "PTA" & contri_boots_scan2$Sample == Sample]) / 
                           sum(contri_boots_unfiltered$mean[contri_boots_unfiltered$sig != "PTA" & contri_boots_unfiltered$Sample == Sample]),
                         Caller = "SCAN2")
  
  PrecRecall <- rbind(PrecRecall, SampleDF_SCAN2)
}

PrecRecall_m <- melt(PrecRecall)

PrecRecall_means <- PrecRecall_m[!is.na(PrecRecall_m$value),] %>% group_by(Caller, variable) %>% summarise(across(value, list(mean = mean, sd = sd)))
PrecRecall_m$Treatment <- factor(PrecRecall_m$Treatment, levels = c("VHC\nNone","MAN\nLow", "MAN\nModerate", "ENU\nLow","ENU\nModerate", "ENU\nHigh"))
S7E <- ggplot(PrecRecall_m, aes(x =variable , y = value, fill = Caller)) + geom_bar(stat = "identity", position = position_dodge(), col = "black", linewidth = 0.4) + facet_grid(.~Treatment, space = "free", scales = 'free') + 
   theme_bw(base_size = 6) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.key.size = unit(0.3, 'cm')) + scale_fill_brewer(palette = "Set1") + scale_y_continuous(expand = c(0,0), limits = c(0,1.05))
 
pdf(file = paste(Output_dir, "FigureS7E_PrecRecall.pdf",sep = ""), width = 3.6, height = 1.7, onefile = F)

  plot_adjusted_strips(S7E, strip_colors = c(Treatment_colors[c(1,2,2,3,3,3)]))
dev.off()


sum(contri_boots_filtered$mean[contri_boots_filtered$sig != "PTA" & contri_boots_scan2$Sample %in% c("ENU-3", "VHC-0")]) / 
  sum(contri_boots_unfiltered$mean[contri_boots_unfiltered$sig != "PTA" & contri_boots_scan2$Sample %in% c("ENU-3", "VHC-0")])

sum(contri_boots_scan2$mean[contri_boots_scan2$sig  != "PTA" & contri_boots_scan2$Sample %in% c("ENU-3", "VHC-0")]) / 
  sum(contri_boots_unfiltered$mean[contri_boots_unfiltered$sig  != "PTA" & contri_boots_scan2$Sample %in% c("ENU-3", "VHC-0")])

# Precision (how many variants are PTA artefacts?)
sum(contri_boots_filtered$mean[contri_boots_filtered$sig == "PTA" & contri_boots_scan2$Sample %in% c("ENU-3", "VHC-0")]) / 
  sum(contri_boots_filtered$mean[contri_boots_filtered$sig != "PTA" & contri_boots_scan2$Sample %in% c("ENU-3", "VHC-0")])

sum(contri_boots_scan2$mean[contri_boots_scan2$sig  == "PTA" & contri_boots_scan2$Sample %in% c("ENU-3", "VHC-0")]) / 
  sum(contri_boots_unfiltered$mean[contri_boots_unfiltered$sig  != "PTA" & contri_boots_scan2$Sample %in% c("ENU-3", "VHC-0")])

## difference in detected SBS5/ENU-mutations between PTATO - SCAN2
1 -  sum(contri_boots_scan2$mean[contri_boots_scan2$sig  %in% c("SBS5", "ENU") & contri_boots_scan2$Sample %in% c("ENU-3", "VHC-0")]) /
   sum(contri_boots_filtered$mean[contri_boots_filtered$sig  %in% c("SBS5", "ENU") & contri_boots_filtered$Sample %in% c("ENU-3", "VHC-0")])
# SCAN2 detects 35% less SBS/ENU mutations


### Indels

# Determine file names
indels_unfilt_vcfs_fnames = list.files(file.path(paste(Input_dir, "GonzalezPena2021_CordBlood/PTATO/", sep= ""), "indels", "Cord"), pattern = ".vcf.gz$", full.names = TRUE)
indels_unfilt_sample_names = get_sample_names(indels_unfilt_vcfs_fnames, str_remove = ".indels.ptato.*")
names(indels_unfilt_vcfs_fnames) <- indels_unfilt_sample_names
# indels_filt_vcfs_fnames = list.files(file.path(Input_dir, "snvs", "Cord"), pattern = "filtered.vcf.gz$", full.names = TRUE, recursive = TRUE)
# snv_filt_sample_names = get_sample_names(snv_filt_vcfs_fnames, str_remove = ".snvs.ptato.*")
# names(snv_filt_vcfs_fnames) <- snv_filt_sample_names

indels_raw <- list()
for(sample in indels_unfilt_sample_names){
  print(sample)
  
  indels_raw[[sample]] <- readIndels(vcf = indels_unfilt_vcfs_fnames[sample], chromosomes = c(1:22), 
                                     VAF_threshold = 0.25)
                                       
}

indel_recurrency_vcf <- readVcf(paste(Input_dir, "Indels/excludelist_v1.vcf", sep = ""), "hg38")
indel_recurrency_grl <- rowRanges(indel_recurrency_vcf)
indel_recurrency_grl <- indel_recurrency_grl[seqnames(indel_recurrency_grl) %in% c(1:22, "X", "Y"),]
seqlevels(indel_recurrency_grl) <- seqlevels(indel_recurrency_grl)[1:24] # The decoy chromosomes give errors later on
seqlevels(indel_recurrency_grl) <- paste("chr", seqlevels(indel_recurrency_grl), sep = "")
indel_recurrency_grl <- get_mut_type(indel_recurrency_grl, type = "indel") # there appears to be one SNVs according to MutationalPatterns

### Characterize the indels in the blacklist
#indel_recurrency_grl <- get_indel_context(indel_recurrency_grl, ref_genome)

indels_raw <- get_mut_type(indels_raw, type = "indel")
indels_raw <- get_indel_context(indels_raw, BSgenome.Hsapiens.UCSC.hg38)

indels_flagged <- filter_indels(indel_grl = indels_raw, indel_recurrency_grl = indel_recurrency_grl, remove_homopolymer = TRUE, mode = "flag")
indels_filtered <- filter_indels(indel_grl = indels_raw, indel_recurrency_grl = indel_recurrency_grl, remove_homopolymer = TRUE)

indel_counts_raw <- count_indel_contexts(indels_raw)
indel_counts_filtered <- count_indel_contexts(indels_filtered)

indel_colors <- c("#FDBE6F", "#FF8001", "#B0DD8B", "#36A12E", "#FDCAB5", 
                  "#FC8A6A", "#F14432", "#BC141A", "#D0E1F2", "#94C4DF", 
                  "#4A98C9", "#1764AB", "#E2E2EF", "#B6B6D8", "#8683BD", 
                  "#61409B")

# unfiltered indels
indel_counts_long <- indel_counts_raw %>% as.data.frame() %>% tibble::rownames_to_column("muttype_total") %>% 
  tidyr::separate(muttype_total, c("muttype", "muttype_sub"), 
                  sep = "_(?=[0-9])") %>% dplyr::mutate(muttype = factor(muttype, 
                                                                         levels = unique(muttype)))
indel_counts_main <- indel_counts_long %>% dplyr::select(-muttype_sub) %>% 
  dplyr::group_by(muttype) %>% dplyr::summarise_all(list(~sum(.))) %>% 
  tidyr::gather(key = "sample", value = "count", -.data$muttype) %>% 
  dplyr::mutate(sample = factor(sample, levels = unique(sample)))
nr_indels <- indel_counts_main %>% dplyr::group_by(sample) %>% dplyr::summarise(nr_muts = sum(count))
facet_labs_y_indels <- stringr::str_c(nr_indels$sample, " (n = ", 
                                      nr_indels$nr_muts, ")")
names(facet_labs_y_indels) <- nr_indels$sample

indel_counts_main$label <- ifelse(indel_counts_main$count > 200, indel_counts_main$count, "")

levels(indel_counts_main$muttype) <- gsub(levels(indel_counts_main$muttype), pattern = "deletion", replacement = "del")
levels(indel_counts_main$muttype) <- gsub(levels(indel_counts_main$muttype), pattern = "insertion", replacement = "ins")
levels(indel_counts_main$muttype) <- gsub(levels(indel_counts_main$muttype), pattern = "microhomology", replacement = "MH")


indel_counts_main$Group <- str_remove(indel_counts_main$sample, "-[0-9]$")
indel_counts_main$Treatment <- str_remove(indel_counts_main$Group, "-[0-9]$")
indel_counts_main$Dose <- str_remove(indel_counts_main$Group, ".*-")
indel_counts_main$Dose_label <- ifelse(indel_counts_main$Dose == 0, "None", "Low" )
indel_counts_main$Dose_label[indel_counts_main$Dose == 3 & indel_counts_main$Treatment == "MAN"] <- "Moderate"
indel_counts_main$Dose_label[indel_counts_main$Dose == 2 & indel_counts_main$Treatment == "ENU"] <- "Moderate"

indel_counts_main$Dose_label[indel_counts_main$Dose == 3 & indel_counts_main$Treatment == "ENU"] <- "High"

mean_indelRaw_perGroup <- indel_counts_main %>% dplyr::group_by(Group, Treatment, Dose, Dose_label, muttype) %>% dplyr::summarise(mean = round(mean(count),0))
mean_indelRaw_perGroup$Treatment <- factor(mean_indelRaw_perGroup$Treatment, levels = c("VHC", "MAN", "ENU"))
mean_indelRaw_perGroup$Dose_label <- factor(mean_indelRaw_perGroup$Dose_label, levels = c("None", "Low", "Moderate", "High"))


indel_raw_counts_main_plot <- ggplot(mean_indelRaw_perGroup, aes(x = Dose_label, y = mean, fill = muttype)) + 
  geom_bar(stat = "identity",  col = "black", width = 0.75, linewidth = 0.2) + 
  geom_text(aes(label = round(stat(y), 0), group = sample),
           stat = 'summary', fun = sum, vjust = -1, size = 1.6) +
  facet_grid(.~Treatment, space = "free", scales = "free", switch = "x") +
  labs(x = "Treatment", y = "Indels (#)", fill = "Indel type") + 
  scale_fill_manual(values = indel_colors) + 
  scale_y_continuous(expand = expansion(mult = c(0,0.1))) +   
  theme_classic(base_size = 6) + 
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), 
        strip.placement = "outside",
        legend.key.size = unit(0.2, 'cm'), 
        legend.margin = margin(t = 5, b = 5, l= 0, r =0, unit = "pt")) +
  guides(fill=guide_legend(ncol=2))

pdf(file = paste(Output_dir, "FigureS7E_Indel_raw.pdf", sep = ""), width = 3.5, height = 2, onefile = F)
plot_adjusted_strips(indel_raw_counts_main_plot, strip_colors = c(Treatment_colors[c(1,2,3)]
))
dev.off()


## Filtered
indel_counts_long <- indel_counts_filtered %>% as.data.frame() %>% tibble::rownames_to_column("muttype_total") %>% 
  tidyr::separate(muttype_total, c("muttype", "muttype_sub"), 
                  sep = "_(?=[0-9])") %>% dplyr::mutate(muttype = factor(muttype, 
                                                                         levels = unique(muttype)))
indel_counts_main <- indel_counts_long %>% dplyr::select(-muttype_sub) %>% 
  dplyr::group_by(muttype) %>% dplyr::summarise_all(list(~sum(.))) %>% 
  tidyr::gather(key = "sample", value = "count", -.data$muttype) %>% 
  dplyr::mutate(sample = factor(sample, levels = unique(sample)))
nr_indels <- indel_counts_main %>% dplyr::group_by(sample) %>% dplyr::summarise(nr_muts = sum(count))
facet_labs_y_indels <- stringr::str_c(nr_indels$sample, " (n = ", 
                                      nr_indels$nr_muts, ")")
names(facet_labs_y_indels) <- nr_indels$sample

indel_counts_main$label <- ifelse(indel_counts_main$count > 200, indel_counts_main$count, "")

levels(indel_counts_main$muttype) <- gsub(levels(indel_counts_main$muttype), pattern = "deletion", replacement = "del")
levels(indel_counts_main$muttype) <- gsub(levels(indel_counts_main$muttype), pattern = "insertion", replacement = "ins")
levels(indel_counts_main$muttype) <- gsub(levels(indel_counts_main$muttype), pattern = "microhomology", replacement = "MH")


indel_counts_main$Group <- str_remove(indel_counts_main$sample, "-[0-9]$")
indel_counts_main$Treatment <- str_remove(indel_counts_main$Group, "-[0-9]$")
indel_counts_main$Dose <- str_remove(indel_counts_main$Group, ".*-")
indel_counts_main$Dose_label <- ifelse(indel_counts_main$Dose == 0, "None", "Low" )
indel_counts_main$Dose_label[indel_counts_main$Dose == 3 & indel_counts_main$Treatment == "MAN"] <- "Moderate"
indel_counts_main$Dose_label[indel_counts_main$Dose == 2 & indel_counts_main$Treatment == "ENU"] <- "Moderate"

indel_counts_main$Dose_label[indel_counts_main$Dose == 3 & indel_counts_main$Treatment == "ENU"] <- "High"

mean_indel_perGroup <- indel_counts_main %>% dplyr::group_by(Group, Treatment, Dose, Dose_label, muttype) %>% dplyr::summarise(mean = round(mean(count),0))
mean_indel_perGroup$Treatment <- factor(mean_indel_perGroup$Treatment, levels = c("VHC", "MAN", "ENU"))
mean_indel_perGroup$Dose_label <- factor(mean_indel_perGroup$Dose_label, levels = c("None", "Low", "Moderate", "High"))
indel_counts_main_plot <- ggplot(mean_indel_perGroup, aes(x = Dose_label, y = mean, fill = muttype)) + 
  geom_bar(stat = "identity",  col = "black", width = 0.75, linewidth = 0.2) + 
  #geom_text( size = 1.6, position = position_stack(vjust = 0.5)) +
  geom_text(aes(label = round(stat(y), 0), group = sample),
            stat = 'summary', fun = sum, vjust = -1, size = 1.6) +
  facet_grid(.~Treatment, space = "free", scales = "free", switch = "x") +
  labs(x = "Treatment", y = "Indels (#)", fill = "Indel type") + 
  scale_fill_manual(values = indel_colors) + 
  scale_y_continuous(expand = expansion(mult = c(0,0.1))) +   
  theme_classic(base_size = 6) + 
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), 
        strip.placement = "outside",
        legend.key.size = unit(0.2, 'cm'), 
        legend.margin = margin(t = 5, b = 5, l= 0, r =0, unit = "pt")) +
  guides(fill=guide_legend(ncol=2))

pdf(file = paste(Output_dir, "FigureS7F_Indel_filtered.pdf", sep = ""), width = 3.5, height = 2, onefile = F)
plot_adjusted_strips(indel_counts_main_plot, strip_colors = c(Treatment_colors[c(1,2,3)]
))
dev.off()

