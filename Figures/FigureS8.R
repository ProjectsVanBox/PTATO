# Supplemental Figure S8

library(VariantAnnotation)
library(MutationalPatterns)
library(BSgenome)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

# Paths
Project_dir <- "/path/to/MendeleyData/directory/"
Metadata_file <-  paste(Project_dir, "Table_S1.txt", sep = "")
Input_dir <- paste(Project_dir, "STE0072/PTATO/snvs/STE0072/", sep = "")
Output_dir <-  "/path/to/output/directory/"
Script_dir <-  "/path/to/GithubScripts/directory/"

if(dir.exists(Output_dir) == F){
  dir.create(Output_dir)
}

source(paste(Script_dir, "GeneralFunctions.R", sep = ""))

# Read the metadata file
Metadata <- read.delim(Metadata_file, header = T)
Metadata_Organoids <- Metadata[Metadata$Individual == "STE0072",]

# Organoid mutation spectra
Organoid_MutSpectra <- read.table(paste(Project_dir, "Resources/Blokzijl_2016.txt",sep = ""),header=T)
STE0072 <- Organoid_MutSpectra[,which(grepl("STE0072", colnames(Organoid_MutSpectra)))]

# PTA artefacts signature
pta_signatures <- read.table(paste(Project_dir, "Resources/PTA_Artefact_Signature.txt", sep = ""),header=T)

# Read the PTATO files
MinimalVAF <- 0.2
SBSs_raw <- list()
for(Sample in Metadata_Organoids$Sample[Metadata_Organoids$Germline == F]){
  print(Sample)
  
  ID <- Metadata$Label[Metadata$Sample == Sample]
  
  vcf_unfiltered <- list.files(paste(Input_dir, sep = ""),
                               pattern = paste(Sample, ".snvs.ptato.vcf.gz$", sep = ""), full.names = T, recursive = T)
  print(vcf_unfiltered)
  vcf_filtered <- list.files(paste(Input_dir, sep = ""),
                             pattern = paste(Sample, ".snvs.ptato.filtered.vcf.gz$", sep = ""), full.names = T, recursive = T)
  print(vcf_filtered)
  SBSs_raw[[ID]] <- readPTATOvcf(vcf_unfiltered = vcf_unfiltered, 
                                     vcf_filtered = vcf_filtered, 
                                     VAF_threshold = MinimalVAF)
}

# Select the right mutations
SBSs_autosomal <- lapply(SBSs_raw, function(x) x[seqnames(x) %in% paste("chr", c(1:22), sep = "")])
SBSs_before <- lapply(SBSs_autosomal, function(x) x[which(x$FILTER != "FAIL_VAF"),])
SBSs_after<- lapply(SBSs_autosomal, function(x) x[which(x$FILTER ==  "PASS"),])

SBSs_removed <- lapply(SBSs_autosomal, function(x) x[which(x$FILTER ==  "FAIL"),])

# Make mutation matrices
mut_mat_before <- mut_matrix(vcf_list = SBSs_before, ref_genome = ref_genome)
mut_mat_after <- mut_matrix(vcf_list = SBSs_after, ref_genome = ref_genome)
mut_mat_removed <- mut_matrix(vcf_list = SBSs_removed, ref_genome = ref_genome)

pdf(file = paste(Output_dir, "FigureS8A_ORG_96profile_before.pdf", sep = ""), width = 3.5, height = 3, onefile = F)
plot_96_profile_small(mut_mat_before, ymax = 0.1, condensed = T)
dev.off()

pdf(file = paste(Output_dir, "FigureS8B_ORG_96profile_after.pdf", sep = ""), width = 3.5, height = 3, onefile = F)
plot_96_profile_small(mut_mat_after, ymax = 0.1, condensed = T)
dev.off()

plot_96_profile(mut_mat_removed)

# Signature refit
invitro_mutmat <- STE0072[,grep("invitro", colnames(STE0072))]

# Change absolute to relative contribution
invitro_mutmat_relative <- sweep(invitro_mutmat,2,colSums(invitro_mutmat),`/`)

# Calculate the mean contribution for each mutation channel for the three organoid samples
invitro_mutmat_relative$Invitro <- rowMeans(invitro_mutmat_relative)

# Merge the organoid and PTA artefact signature
signatures <- cbind(PTA = pta_signatures[,1], In_vitro = invitro_mutmat_relative$Invitro)
rownames(signatures) <- rownames(pta_signatures)

colnames(signatures) <- c("PTA", "Organoid\n(In vitro)")
pdf(file = paste(Output_dir, "FigureS8C_ORG_Signatures.pdf", sep = ""), width = 3.5, height = 1.5, onefile = F)
plot_96_profile_small(signatures, ymax = 0.1, condensed = T)
dev.off()

### Perform bootstrapped signature refitting
## Before PTATO
contri_boots_before <- fit_to_signatures_bootstrapped(mut_mat_before,
                                                      signatures,
                                                      n_boots = 100,
                                                      method = "strict", max_delta = 0.005)
# Change to long format for plotting
contri_before_tb <- contri_boots_before %>% as.data.frame() %>% tibble::rownames_to_column("exp") %>% 
  tidyr::gather(key = "sig", value = "contri", -exp) %>% 
  dplyr::mutate(Sample = gsub("_[^_]+$", "", exp), Sample = factor(Sample, 
                                                                   levels = unique(Sample)), sig = factor(sig, levels = unique(sig)))
# Calculate the mean contribution for each signature
contri_before_tb2 <- contri_before_tb %>% dplyr::group_by(Sample, sig) %>% 
  dplyr::summarise(mean = mean(contri)) %>% 
  dplyr::ungroup()
contri_before_tb2$Type <- "Before PTATO"

## After PTATO
contri_boots_after <- fit_to_signatures_bootstrapped(mut_mat_after,
                                                      signatures,
                                                      n_boots = 100,
                                                      method = "strict", max_delta = 0.005)
# Change to long format for plotting
contri_after_tb <- contri_boots_after %>% as.data.frame() %>% tibble::rownames_to_column("exp") %>% 
  tidyr::gather(key = "sig", value = "contri", -exp) %>% 
  dplyr::mutate(Sample = gsub("_[^_]+$", "", exp), Sample = factor(Sample, 
                                                                   levels = unique(Sample)), sig = factor(sig, levels = unique(sig)))
# Calculate the mean contribution for each signature
contri_after_tb2 <- contri_after_tb %>% dplyr::group_by(Sample, sig) %>% 
  dplyr::summarise(mean = mean(contri)) %>% 
  dplyr::ungroup()
contri_after_tb2$Type <- "After PTATO"

## Removed by PTATO
contri_boots_removed <- fit_to_signatures_bootstrapped(mut_mat_removed,
                                                     signatures,
                                                     n_boots = 100,
                                                     method = "strict", max_delta = 0.005)

# Change to long format for plotting
contri_removed_tb <- contri_boots_removed %>% as.data.frame() %>% tibble::rownames_to_column("exp") %>% 
  tidyr::gather(key = "sig", value = "contri", -exp) %>% 
  dplyr::mutate(Sample = gsub("_[^_]+$", "", exp), Sample = factor(Sample, 
                                                                   levels = unique(Sample)), sig = factor(sig, levels = unique(sig)))
# Calculate the mean contribution for each signature
contri_removed_tb2 <- contri_removed_tb %>% dplyr::group_by(Sample, sig) %>% 
  dplyr::summarise(mean = mean(contri)) %>% 
  dplyr::ungroup()
contri_removed_tb2$Type <- "Removed by PTATO"


## Merge refit data
Refits_merged <- rbind(contri_before_tb2, contri_after_tb2, contri_removed_tb2)
Refits_merged$Type <- factor(Refits_merged$Type, levels = c("Before PTATO","Removed by PTATO","After PTATO"))

# Rename signature
Refits_merged$sig <- factor(Refits_merged$sig, levels = c("PTA", "Organoid\n(In vitro)"))

ggplot(Refits_merged, aes(x = Sample, y = mean, fill = sig)) + geom_bar(stat = "identity", width = 0.9) + 
  facet_grid(.~Type) + 
  theme_classic(base_size = 6) +
  scale_fill_manual(values = c( "#D55E00","#009E73")) + scale_y_continuous(expand = c(0,0)) +
  labs(y = "Single base\nsubstitutions (#)", fill = "Mutational profile") +
  theme(legend.key.size = unit(0.4, 'cm'),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste(Output_dir, "FigureS8D_ORG_Refit.pdf", sep = ""), width = 9, height = 4, units = "cm")

