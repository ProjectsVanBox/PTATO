### Figures 2, S5 and S6

library(VariantAnnotation)
library(MutationalPatterns)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
library(reshape2)
library(ggh4x) # For nested facets (indel profiles)

# Input dir should contain the raw PTATO/SMuRF/CallableLoci output, eg: /[INPUT_DIR]/[CELL_LINE]/PTATO/
Input_dir <- "/path/to/MendeleyData/directory/"
Script_dir <- "/path/to/GithubScripts/directory/"
# Output figures and data will be written in the Output_Dir
Output_dir <- "/path/to/output/directory/"

Output_dir_Figure2 <- paste(Output_dir, "Figure2/", sep = "")

Output_dir_FigureS5 <- paste(Output_dir, "FigureS5/", sep = "")
Output_dir_FigureS6 <- paste(Output_dir, "FigureS6/", sep = "")

if(dir.exists(Output_dir_Figure2) == F){dir.create(Output_dir_Figure2)}
if(dir.exists(Output_dir_FigureS5) == F){dir.create(Output_dir_FigureS5)}
if(dir.exists(Output_dir_FigureS6) == F){dir.create(Output_dir_FigureS6)}

source(paste(Script_dir, "GeneralFunctions.R", sep = ""))

# The metadata files contains the sample names, labels, and days in culture
Metadata_file <- paste(Input_dir, "Table_S2.txt", sep = "")
if(file.exists(Metadata_file) == T){
  print(Metadata_file)
  Metadata <- read.delim(Metadata_file, header = T)
}

# Color settings
AHH1_colors <- c("#8DA0CB", "#66C2A5","#FC8D62")
names(AHH1_colors) <- c("WT", "FANCCKO", "MSH2KO")

sSNVs_raw <- list()
sIndels_PTA_raw <- list()
sIndels_Clones <- list()

sSNVs_PASS <- list()
sSNVs_FAIL <- list()

sSNVs_absentClone_raw <- list()
sSNVs_absentClone <- list()
sIndels_absentClone_raw <- list()
sIndels_absentClone <- list()

# These lists will contain the variants that are not present in the preceding (sub)clones (they may be present in other samples derived from this sample)
sSNVs_unique_raw <- list()
sSNVs_unique <- list()
sIndels_unique_raw <- list()
sIndels_unique <- list()

sSNVs_SCAN2 <- list()
sSNVs_SCAN2_unique <- list()


### 
# Read the CallableLoci files (takes a while)
CallableLoci <- list()
for(Sample in Metadata$Sample[ Metadata$Type == "PTA"]){
  print(Sample)
  
  CallableLoci_file <- list.files(paste(Input_dir, Metadata$Cell_line[Metadata$Sample == Sample], "/CallableLoci/", sep = ""), pattern = paste(Sample, "_CallableLoci.bed", sep = "" ), full.names = T)
  print(CallableLoci_file)
  print(paste("# Reading: ", CallableLoci_file, sep = ""))
  CallableLoci_sample <- read.delim(CallableLoci_file, header = F)
  Callable <- CallableLoci_sample[CallableLoci_sample[,4] == "CALLABLE",]
  CallableLoci_g <- GRanges(seqnames = Callable[,1], IRanges(start = Callable[,2], end = Callable[,3]))
  CallableLoci[[Sample]] <- CallableLoci_g
}

# # Read the fraction of the genome that is callable:
# Metadata$CallableLoci <- 0.95
# for(Sample in Metadata$Sample[Metadata$Cell_line != "PMCAHH1"]){
#   print(Sample)
#   callable_file <- list.files(paste(Input_dir, Metadata$Cell_line[Metadata$Sample == Sample], "/PTATO/snvs_callable/", sep = ""), pattern = paste(Sample, "_CallableLoci.txt$", sep = ""), full.names = T)
#   if(length(callable_file) == 1){
#     print(callable_file)
#     callable_sample <- read.delim(callable_file)
#     Metadata$CallableLoci[Metadata$Sample == Sample] <- callable_sample$Freq[callable_sample$State == "CALLABLE"]
#   }
# }

# Read the PTATO VCFs for PTA samples and the SMuRF VCFs for clones/bulk
# Exclude PMCAHH1, which is the bulk control sample

MinimalVAF <- 0.2 # all variants show have a VAF higher than MinimalVAF

for(Sample in Metadata$Sample[Metadata$Cell_line != "PMCAHH1"]){
  print(Sample)
  if(Metadata$Type[Metadata$Sample == Sample] == "PTA"){
    snv_vcf_file <-  list.files(paste(Input_dir, Metadata$Cell_line[Metadata$Sample == Sample], "/PTATO/snvs_callable/", sep = ""), pattern = paste(Sample, ".snvs.ptato.callable.vcf.gz$", sep = ""), full.names = T)
    print(snv_vcf_file)
    snv_vcf_filtered_file <-  list.files(paste(Input_dir, Metadata$Cell_line[Metadata$Sample == Sample], "/PTATO/snvs/", sep = ""), pattern = paste(Sample, ".snvs.ptato.filtered.vcf.gz$", sep = ""), full.names = T, recursive = T)
    print(snv_vcf_filtered_file)
    
    indel_vcf_file <-  list.files(paste(Input_dir, Metadata$Cell_line[Metadata$Sample == Sample], "/PTATO/indels/", sep = ""), pattern = paste(Sample, ".indels.ptato.vcf.gz$", sep = ""), full.names = T, recursive = T)
    print(indel_vcf_file)
    
    if(length(snv_vcf_file) == 1 & length(snv_vcf_filtered_file) == 1 & length(indel_vcf_file) == 1){
      # The filtered VCF is necessary to retrieve the PTAprobsCutoffs
      sSNVs_raw[[Sample]] <- readPTATOvcf(vcf_unfiltered = snv_vcf_file, 
                                            vcf_filtered = snv_vcf_filtered_file, 
                                            VAF_threshold = MinimalVAF, 
                                          cutoff = 2, chromosomes = c(1:22))
      # Some variants have a SNP-ID as the variant name, instead of the genomic coordinates
      names(sSNVs_raw[[Sample]]) <- paste(gsub(pattern = "chr", replacement = "", as.vector(seqnames(sSNVs_raw[[Sample]]))), ":", start(sSNVs_raw[[Sample]]), "_",  as.vector(sSNVs_raw[[Sample]]$REF), "/",as.vector(unlist(sSNVs_raw[[Sample]]$ALT)), sep = "")
      

      sIndels_PTA_raw[[Sample]] <- readIndels(indel_vcf_file, chromosomes = c(1:22), VAF_threshold = MinimalVAF, rename_indels = T)
      

    }
  } else if (Metadata$Type[Metadata$Sample == Sample] == "Subclone"){
    
    vcf_file <- list.files(paste(Input_dir, Metadata$Cell_line[Metadata$Sample == Sample], "/PTATO/intermediate/short_variants/somatic_vcfs/",Metadata$Cell_line[Metadata$Sample == Sample],"/", sep = ""), pattern = paste(Sample, ".SMuRF.filtered.vcf.gz$", sep = ""), full.names = T)
    print(vcf_file)
    vcf <- readVcf(vcf_file, genome = "hg38")
    vcf <- vcf[as.vector(seqnames(vcf)) %in% c(1:22),]
    snv_vcf <- vcf[isSNV(vcf),]
    
    GR <- VCF_to_GR(snv_vcf)
    GR$VAF <- as.numeric(geno(snv_vcf)$VAF)
    GR$PTAprob <- NA
    GR$PTAprobCutoff <- NA 
    GR$PTAprobCutoff_Walker <- NA
    GR$PTAprobCutoff_Cossim <- NA
    names(GR) <- paste(gsub(pattern = "chr", replacement = "", x = as.vector(seqnames(GR))), ":", start(GR), "_",  
                        as.vector(GR$REF), "/",as.vector(unlist(GR$ALT)), sep = "")
    
    sSNVs_raw[[Sample]] <- GR[GR$VAF > MinimalVAF]
    
    # Indels
    sIndels_Clones[[Sample]] <-  readIndels(vcf_file, chromosomes = c(1:22), VAF_threshold = MinimalVAF, rename_indels = T)
  }
}

# Remove all variants that are UNCALLABLE and/or with a PTAprobs < PTAprobsCutoff
sSNVs_PASS <- lapply(sSNVs_raw, function(x) x[which(x$FILTER =="PASS"),])
sSNVs_FAIL <- lapply(sSNVs_raw, function(x) x[which(x$FILTER =="FAIL"),])

### Read the exclusion list
### Reading with Mutationalpatterns doesn't work well (it thinks many indels at same position are duplicates, even though they have different ALT)
#blacklist_grl <- readVcf("~/hpc/pmc_vanboxtel/personal/smiddelkamp/PTA_Indel/VCF4/PTA_Indel_Blacklist_normNoGT.vcf.gz", "hg38")
exclusion_vcf <- readVcf(paste(Input_dir, "Indels/excludelist_v1.vcf", sep = ""), "hg38")
indel_recurrency_grl <- rowRanges(exclusion_vcf)
indel_recurrency_grl <- indel_recurrency_grl[seqnames(indel_recurrency_grl) %in% c(1:22, "X", "Y"),]
seqlevels(indel_recurrency_grl) <- seqlevels(indel_recurrency_grl)[1:24] # The decoy chromosomes give errors later on
seqlevels(indel_recurrency_grl) <- paste("chr", seqlevels(indel_recurrency_grl), sep = "")
indel_recurrency_grl <- get_mut_type(indel_recurrency_grl, type = "indel") # there appear to be SNVs according to MutationalPatterns


RecurrencyList_Info <- read.table(paste(Input_dir, "Indels/Indel_RecurrencyList_Info.txt", sep = ""), sep = ";")
RecurrencyList_Info$V3 <- gsub(RecurrencyList_Info$V3, pattern = "IC=", replacement = "")

exclusion3ind <- exclusion_vcf[which(RecurrencyList_Info$V3 > 2),]
exclusion3ind_grl <- rowRanges(exclusion3ind)
exclusion3ind_grl <- exclusion3ind_grl[seqnames(exclusion3ind_grl) %in% c(1:22, "X", "Y"),]
seqlevels(exclusion3ind_grl) <- seqlevels(exclusion3ind_grl)[1:24] # The decoy chromosomes give errors later on
seqlevels(exclusion3ind_grl) <- paste("chr", seqlevels(exclusion3ind_grl), sep = "")
exclusion3ind_grl <- get_mut_type(exclusion3ind_grl, type = "indel") # there appears to be one SNVs according to MutationalPatterns


sIndels_PTA_raw <- get_mut_type(sIndels_PTA_raw, type = "indel")
sIndels_PTA_raw <- get_indel_context(sIndels_PTA_raw, ref_genome)

sIndels_PTA_flag <- filter_indels(indel_grl = sIndels_PTA_raw, indel_recurrency_grl = indel_recurrency_grl, remove_homopolymer = TRUE, mode = "flag")
sIndels_PTA_filter <- filter_indels(indel_grl = sIndels_PTA_raw, indel_recurrency_grl = indel_recurrency_grl, remove_homopolymer = TRUE)

sIndels_PTA_flag <- filter_indels(indel_grl = sIndels_PTA_raw, indel_recurrency_grl = exclusion3ind_grl, remove_homopolymer = TRUE, mode = "flag")
sIndels_PTA_filter <- filter_indels(indel_grl = sIndels_PTA_raw, indel_recurrency_grl = exclusion3ind_grl, remove_homopolymer = TRUE)


sIndels_raw <- c(sIndels_PTA_flag, sIndels_Clones)
sIndels_merged <- c(sIndels_PTA_filter, sIndels_Clones)

### Read the SCAN2 output
SCAN2_AHH1_Bulk <- read.delim(paste(Input_dir,"PMCAHH1-WT/SCAN2/rescued_muts.txt", sep = ""))
SCAN2_FANCC_MSH2_Bulk <- read.delim(paste(Input_dir,"PMCAHH1-FANCCKO/SCAN2/rescued_muts.txt", sep = ""))

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

sSNVs_SCAN2 <- c(SCAN2_to_GRanges(rescued_muts = SCAN2_AHH1_Bulk, type = "snv"), SCAN2_to_GRanges(rescued_muts = SCAN2_FANCC_MSH2_Bulk, type = "snv"))

### Read the multisample VCFs and determine which variants are present in which samples (variant_tables)
multisample_VCFs <- list()
variant_tables <- list()

for(Cell_line in unique(Metadata$Cell_line[Metadata$Cell_line != "PMCAHH1"])){
  print(Cell_line)
  multisample_vcf_file <- list.files(paste(Input_dir, Cell_line, "/PTATO/intermediate/short_variants/SMuRF/", sep = ""), 
                                     pattern = ".SMuRF.filtered.vcf$", full.names = T, recursive = T)
  print(multisample_vcf_file)
  multisample_vcf <- readVcf(multisample_vcf_file, genome = "hg38")
  
  # Only select variants on the autosomes:
  multisample_vcf <- multisample_vcf[as.vector(seqnames(multisample_vcf) %in% c(1:22))]
  
  names(multisample_vcf) <- paste(as.vector(seqnames(rowRanges(multisample_vcf))), ":", start(rowRanges(multisample_vcf)), "_",  as.vector(rowRanges(multisample_vcf)$REF), "/",as.vector(unlist(rowRanges(multisample_vcf)$ALT)), sep = "")
  
  multisample_VCFs[[Cell_line]] <- multisample_vcf     
  
  # Count in which cell a variant is present (=1) or absent (=0)
  variant_table <- matrix(nrow=length(rowRanges(multisample_vcf)),ncol=ncol(multisample_vcf))
  colnames(variant_table) <- colnames(multisample_vcf)
  
  for ( sample in colnames(variant_table)) {
    tmp_table <- unlist(lapply(info(multisample_vcf)$CLONAL_SAMPLE_NAMES, function(x) match(sample, x)))
    tmp_table2 <- unlist(lapply(info(multisample_vcf)$FAIL_QC_SAMPLE_NAMES, function(x) match(sample, x)))
    
    tmp_table[is.na(tmp_table)] <- 0
    tmp_table[tmp_table > 0] <- 1
    # Sample that FAIL_QC get a 1 in tmp_table2. So everything that's is a 1 (and not na) should be removed
    tmp_table[!is.na(tmp_table2)] <- NA
    variant_table[,sample] <- tmp_table
  }
  
  # #The names in the VCF contain SNP ids, the names in the variant_table do not
  # row.names(variant_table) <- paste(as.vector(seqnames(rowRanges(multisample_vcf))), ":", start(rowRanges(multisample_vcf)), "_",  as.vector(rowRanges(multisample_vcf)$REF), "/",as.vector(unlist(rowRanges(multisample_vcf)$ALT)), sep = "")
  
  row.names(variant_table) <- names(multisample_vcf)
  colSums(variant_table, na.rm = T)
  
  variant_tables[[Cell_line]] <- variant_table
}

# Remove the variants in the clone 
for(Cell_line in Metadata$Cell_line[Metadata$Cell_line != "PMCAHH1"]){
  print(Cell_line)
  # Select the multi-sample variant_table for this cell line
  variant_table <- variant_tables[[Cell_line]] 
  for(Sample in Metadata$Sample[Metadata$Cell_line == Cell_line]){
    print(Sample)
    Clone <- Metadata$Clone[Metadata$Sample == Sample]
    #print(Clone)
    # Determine which variants are present in the clone
    Variants_in_clone <- variant_table[which(variant_table[,Clone] != 0),]

    if(nrow(Variants_in_clone) > 0){
      sSNVs_absentClone_raw[[Sample]] <- sSNVs_raw[[Sample]][-which(names(sSNVs_raw[[Sample]]) %in% row.names(Variants_in_clone))]
      sSNVs_absentClone[[Sample]] <- sSNVs_PASS[[Sample]][-which(names(sSNVs_PASS[[Sample]]) %in% row.names(Variants_in_clone))]
      sIndels_absentClone_raw[[Sample]] <- sIndels_raw[[Sample]][-which(names(sIndels_raw[[Sample]]) %in% row.names(Variants_in_clone))]
      sIndels_absentClone[[Sample]] <- sIndels_merged[[Sample]][-which(names(sIndels_merged[[Sample]]) %in% row.names(Variants_in_clone))]
      
      #SNV_SCAN2_AbsentClone <- SNV_SCAN2[[Sample]][-which(names(SNV_SCAN2[[Sample]]) %in% row.names(Variants_in_clone))]
    } else {
      sSNVs_absentClone_raw[[Sample]] <- sSNVs_raw[[Sample]]
      sSNVs_absentClone[[Sample]] <- sSNVs_PASS[[Sample]]
      sIndels_absentClone_raw[[Sample]] <- sIndels_raw[[Sample]]
      sIndels_absentClone[[Sample]] <- sIndels_merged[[Sample]]
      #SNV_SCAN2_AbsentClone[[Sample]] <- SNV_SCAN2[[Sample]]
    }
    
    # The unique variants contain the variants that arose after the clonal step preceding the sample (eg absent in all the (sub)clones before the sample). They can be shared with daughter samples
    Variants_preceding_sample <- variant_table[,c(1:which(colnames(variant_table) == Sample))]
    if(!is.null(nrow(Variants_preceding_sample))){
      Unique_variants <- Variants_preceding_sample[which(Variants_preceding_sample[,Sample] == 1 & rowSums(Variants_preceding_sample, na.rm = T) == 1),]
    } else {
      Unique_variants <- Variants_preceding_sample[which(Variants_preceding_sample == 1 & sum(Variants_preceding_sample, na.rm = T) == 1)]
      
    }
    if(!is.null(nrow(Variants_preceding_sample))){
      sSNVs_unique_raw[[Sample]] <- sSNVs_raw[[Sample]][which(names(sSNVs_raw[[Sample]]) %in% row.names(Unique_variants))]
      sSNVs_unique[[Sample]] <- sSNVs_PASS[[Sample]][which(names(sSNVs_PASS[[Sample]]) %in% row.names(Unique_variants))]
      sIndels_unique_raw[[Sample]] <- sIndels_raw[[Sample]][which(names(sIndels_raw[[Sample]]) %in% row.names(Unique_variants))]
      sIndels_unique[[Sample]] <- sIndels_merged[[Sample]][which(names(sIndels_merged[[Sample]]) %in% row.names(Unique_variants))]
      
      if(Sample %in% names(sSNVs_SCAN2)){
        sSNVs_SCAN2_unique[[Sample]] <- sSNVs_SCAN2[[Sample]][which(names(sSNVs_SCAN2[[Sample]]) %in% row.names(Unique_variants))]
      }
      
    }
  }
}



# library(pheatmap)
# pheatmap(variant_table, show_rownames = F)
# pheatmap(variant_tables[["PMCAHH1-FANCCKO"]], show_rownames = F)

### SNVs per day since first clonal step
VariantCounts <- data.frame(Sample = names(sSNVs_absentClone), 
                            PTATO = lengths(sSNVs_absentClone), 
                            Raw = lengths(sSNVs_absentClone_raw))
VariantCounts <- merge(VariantCounts, Metadata, by = "Sample")

VariantCounts_m <- melt(VariantCounts[,c("Label", "Days_after_clone", "Gene", "PTATO", "Raw", "Type", "CallableLoci")], 
                        id.vars = c("Label", "Days_after_clone", "Gene", "Type", "CallableLoci"))
VariantCounts_m$variable <- as.character(VariantCounts_m$variable)
VariantCounts_m$Muts_Callable <- VariantCounts_m$value / VariantCounts_m$CallableLoci
# For the subclones we only show the raw (not PTATO filtered) number of sSNVs
VariantCounts_m$variable[VariantCounts_m$Type != "PTA"] <- "Clone"
VariantCounts_m <- VariantCounts_m[!duplicated(paste(VariantCounts_m$Label, VariantCounts_m$variable)),]

# The non-PTA samples will be used to calculate the expected mutation rate and burden, which will be shown as the ablines
VariantCounts_nonPTA <- VariantCounts_m[VariantCounts_m$Type != "PTA" ,]
VariantCounts_nonPTA$Muts_per_day <- VariantCounts_nonPTA$Muts_Callable/VariantCounts_nonPTA$Days_after_clone
# Per gene knockout, calculate the mean, max and min mutation rate over the subclones
VariantCounts_nonPTA_summary <- as.data.frame(VariantCounts_nonPTA %>% group_by(Gene) %>%
                                                summarise(across(Muts_per_day, list(mean = mean, min = min, max = max))))
Mutation_intervals <- data.frame(x = c(0,0,0, rep(max(VariantCounts$Days_after_clone+5), 3)),
                                 ymin = c(0,0,0,  VariantCounts_nonPTA_summary$Muts_per_day_min*max(VariantCounts$Days_after_clone+5 )),
                                 ymax = c(0,0,0, VariantCounts_nonPTA_summary$Muts_per_day_max*max(VariantCounts$Days_after_clone+5)),
                                 Gene = rep(VariantCounts_nonPTA_summary$Gene, 2))

VariantCounts_m <- rbind(VariantCounts_m, data.frame(Label = c("WT-Clone1", "WT-Clone2", "FANCCKO-Clone1", "MSH2KO-Clone1"),
                                                     Days_after_clone = 0,
                                                     Gene = c("WT", "WT", "FANCC", "MSH2"),
                                                     Type = "Clone", CallableLoci = NA, variable = "Clone",  value = 0, Muts_Callable = 0))

ggplot(VariantCounts_m, aes(x = Days_after_clone, y = Muts_Callable, fill = Gene, shape = variable)) + 
  geom_ribbon(data=Mutation_intervals,
              aes(x=x,ymin=ymin, ymax=ymax, fill = Gene),
               alpha = 0.2, inherit.aes = F) +
  geom_abline(data = VariantCounts_nonPTA_summary, aes(slope = Muts_per_day_mean, intercept = 0, col = Gene), linetype = c(1,2,2)) +
  geom_point(size = 1.5 ) +
  coord_cartesian(xlim = c(0, max(VariantCounts$Days_after_clone+5)), ylim = c(0, max(VariantCounts_m$Muts_Callable)+100)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_shape_manual(values = c(22, 23, 21)) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme_classic(base_size = 6) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth= 0.5)) +
  labs(y = "Base substitution accumulation\n(since first clonal step)", x = "Time since first clonal step (Days)", shape = "Filter")

ggsave(paste(Output_dir_Figure2, "Figure2C_PMCAHH1_SBSAccumulation.pdf", sep = ""), width = 7, height = 4, units = "cm")


### SNVs per day since first clonal stepx
IndelCounts <- data.frame(Sample = names(sIndels_absentClone), 
                            PTATO = lengths(sIndels_absentClone), 
                            Raw = lengths(sIndels_absentClone_raw))
IndelCounts <- merge(IndelCounts, Metadata, by = "Sample")

IndelCounts_m <- melt(IndelCounts[,c("Label", "Days_after_clone", "Gene", "PTATO", "Raw", "Type", "CallableLoci")], 
                        id.vars = c("Label", "Days_after_clone", "Gene", "Type", "CallableLoci"))
IndelCounts_m$variable <- as.character(IndelCounts_m$variable)
IndelCounts_m$Indels_Callable <- IndelCounts_m$value / IndelCounts_m$CallableLoci
# For the subclones we only show the raw (not PTATO filtered) number of sSNVs
IndelCounts_m$variable[IndelCounts_m$Type != "PTA"] <- "Clone"
IndelCounts_m <- IndelCounts_m[!duplicated(paste(IndelCounts_m$Label, IndelCounts_m$variable)),]

# The non-PTA samples will be used to calculate the expected mutation rate and burden, which will be shown as the ablines
IndelCounts_m_nonPTA <- IndelCounts_m[IndelCounts_m$Type != "PTA" ,]
IndelCounts_m_nonPTA$Indels_per_day <- IndelCounts_m_nonPTA$Indels_Callable/ IndelCounts_m_nonPTA$Days_after_clone
# Per gene knockout, calculate the mean, max and min mutation rate over the subclones
IndelCounts_nonPTA_summary <- as.data.frame(IndelCounts_m_nonPTA %>% group_by(Gene) %>%
                                                summarise(across(Indels_per_day, list(mean = mean, min = min, max = max))))
Indel_intervals <- data.frame(x = c(0,0,0, rep(max(IndelCounts$Days_after_clone+5), 3)),
                                 ymin = c(0,0,0,  IndelCounts_nonPTA_summary$Indels_per_day_min*max(IndelCounts$Days_after_clone+5 )),
                                 ymax = c(0,0,0, IndelCounts_nonPTA_summary$Indels_per_day_max*max(IndelCounts$Days_after_clone+5)),
                                 Gene = rep(IndelCounts_nonPTA_summary$Gene, 2))

ggplot(IndelCounts_m, aes(x = Days_after_clone, y = Indels_Callable, fill = Gene, shape = variable)) + 
  geom_ribbon(data=Indel_intervals,
              aes(x=x,ymin=ymin, ymax=ymax, fill = Gene),
              alpha = 0.2, inherit.aes = F) +
  geom_abline(data = IndelCounts_nonPTA_summary, aes(slope = Indels_per_day_mean, intercept = 0, col = Gene), linetype = c(1,2,2)) +
  geom_point(size = 1.5 ) +
  coord_cartesian(xlim = c(0, max(IndelCounts$Days_after_clone+5)), ylim = c(0, max(IndelCounts_m$Indels_Callable)+100)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_shape_manual(values = c(22, 23, 21)) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme_classic(base_size = 6) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth= 0.5)) +
  labs(y = "Indel accumulation\n(since first clonal step)", x = "Time since first clonal step (Days)", shape = "Filter")

ggsave(paste(Output_dir_FigureS6, "FigureS6C_PMCAHH1_IndelAccumulation.pdf", sep = ""), width = 7, height = 4, units = "cm")

### Base substitutions observed vs expected
sSNVs_unique_PASS <- lapply(sSNVs_unique_raw, function(x) x[which(x$FILTER =="PASS"),]) 
sSNVs_unique_RAW <- lapply(sSNVs_unique_raw, function(x) x[which(x$FILTER %in% c("PASS","FAIL")),])  # only select the variants with sufficient coverage and PASS_QC (these variants would also be filtered away by standard non-PTATO filtering)
sSNVs_unique_FAIL <- lapply(sSNVs_unique_raw, function(x) x[which(x$FILTER =="FAIL"),]) 

SCAN2_burdens <- read.delim(paste(Input_dir, "PMCAHH1_SCAN2_SNV_Burden.txt", sep = ""))

Unique_counts <- data.frame(Sample = names(sSNVs_unique), 
                            PTATO = lengths(sSNVs_unique_PASS), 
                            BeforePTATO = lengths(sSNVs_unique_RAW),
                            RemovedPTATO = lengths(sSNVs_unique_FAIL)) 
Unique_counts_SCAN2 <- data.frame(Sample =SCAN2_burdens$Sample, SCAN2 = SCAN2_burdens$Extrapolated)
Unique_counts <- merge(Unique_counts, Unique_counts_SCAN2, all.x = T)
Unique_counts <- merge(Unique_counts, Metadata, by = "Sample")

Unique_counts_m <- melt(Unique_counts[,c("Label", "Days_after_sort", "Gene", "PTATO", "BeforePTATO", "RemovedPTATO",  "SCAN2", "Type", "CallableLoci")], 
                        id.vars = c("Label", "Days_after_sort", "Gene", "Type", "CallableLoci"))
Unique_counts_m$Muts_Callable <- Unique_counts_m$value / Unique_counts_m$CallableLoci
# SCAN2 already performs a extrapolation step (so shouldn't be extrapolated also for callable loci)
Unique_counts_m$Muts_Callable[Unique_counts_m$variable == "SCAN2"] <-  Unique_counts_m$value[Unique_counts_m$variable == "SCAN2"] 
Unique_counts_m <- merge(Unique_counts_m, VariantCounts_nonPTA_summary, by = "Gene", all.x = T)

Unique_counts_m$Expected <- Unique_counts_m$Muts_per_day_mean * Unique_counts_m$Days_after_sort
Unique_counts_m$ObsVsExp <- Unique_counts_m$Muts_Callable / Unique_counts_m$Expected

Unique_counts_m$Muts_per_day <- Unique_counts_m$Muts_Callable / Unique_counts_m$Days_after_sort

Unique_counts_m$variable <- as.vector(Unique_counts_m$variable)
Unique_counts_m$variable[Unique_counts_m$variable == "BeforePTATO"] <- "Before\nPTATO"
Unique_counts_m$variable[Unique_counts_m$variable == "RemovedPTATO"] <- "Removed\nby PTATO"

Unique_counts_m$variable <- factor(Unique_counts_m$variable, levels = c("SCAN2", "PTATO","Removed\nby PTATO", "Before\nPTATO"))

Unique_counts_m$Linestart <- 1
Unique_counts_PTA <- Unique_counts_m[Unique_counts_m$Type == "PTA",]
Unique_counts_PTA$Label2 <- gsub(pattern = "-", replacement = "-\n", Unique_counts_PTA$Label)
Unique_counts_PTA$Label2 <- factor(Unique_counts_PTA$Label2, levels = c("MSH2KO-\nPTA1", "FANCCKO-\nPTA1","WT-\nPTA2","WT-\nPTA1" ))

# Print the accuracy of mutation burden prediction
Unique_counts_PTA$Accuracy <- 1-abs(1-Unique_counts_PTA$ObsVsExp)
mean(Unique_counts_PTA$Accuracy[Unique_counts_PTA$variable == "PTATO"])
mean(Unique_counts_PTA$Accuracy[Unique_counts_PTA$variable == "SCAN2"])

Burden <- Unique_counts_PTA %>% group_by(variable) %>%
  summarise(across(c(Accuracy,ObsVsExp), list(mean = mean, sd = sd)))
ggplot(Burden, aes( x= variable, y = ObsVsExp_mean)) + geom_bar(stat = "identity", width = 0.8) + 
  
  geom_errorbar(aes(ymin=ObsVsExp_mean-ObsVsExp_sd, ymax=ObsVsExp_mean+ObsVsExp_sd), width=.2,
                position=position_dodge(.9))  +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.4) +
  geom_text(aes(y = 1.9, label = round(Accuracy_mean,2)), size = 1.8, fontface = 3) +
  coord_flip() + 
  scale_y_continuous(expand = c(0,0.01), limits = c(0,2)) +
  theme_classic(base_size = 6) +
  labs(x = "Unique base substitutions", y= "Base substitution burden\n(Observed vs Expected)")
ggsave(paste(Output_dir_Figure2, "Fig_2c_PMCAHH1_MutationObsVsExp2.pdf", sep = ""), width = 4, height = 3.25, units = "cm")


ggplot(Unique_counts_PTA, 
       aes(x = Label2, y = ObsVsExp, col = variable )) + 
  geom_hline(yintercept = 1, linewidth = 0.2, linetype = 2) + coord_flip() +
  geom_linerange(aes(ymin = Linestart, ymax = ObsVsExp), linetype = 3, linewidth = 0.2) +  
  geom_point(aes(shape  = variable, fill = variable)) +
  scale_y_continuous(limits = c(0,2)) +
  # scale_shape_manual(values = c(23, 21)) +
  # scale_fill_manual(values = c( "#009E73","#D55E00")) +
  # scale_color_manual(values = c("#009E73","#D55E00")) +
   theme_classic(base_size = 6) +
  labs(fill = "Filter", color = "Filter", shape = "Filter", y = "Unique base substitutions\nObserved/Expected", x = "Sample") +
  theme(legend.margin=margin( l = 0, r = 0, unit='cm'),
        panel.border = element_rect(colour = "black", fill=NA, linewidth= 0.5))

#ggsave(paste(Output_dir_Figure2, "Fig_2c_PMCAHH1_MutationObsVsExp.pdf", sep = ""), width = 5.5, height = 3, units = "cm")


### Indels observed vs expected
indel_counts_fail <- count_indel_contexts( lapply(sIndels_unique_raw, function(x) x[which(x$FILTER != "PASS"),]) )

Unique_indel_counts <- data.frame(Sample = names(sIndels_unique), 
                            PTATO = lengths(sIndels_unique), 
                            BeforePTATO = lengths(sIndels_unique_raw),
                            RemovedPTATO = lengths( lapply(sIndels_unique_raw, function(x) x[which(x$FILTER != "PASS"),]))) 
Unique_indel_counts <- merge(Unique_indel_counts, Metadata, by = "Sample")

Unique_indel_counts_m <- melt(Unique_indel_counts[,c("Label", "Days_after_sort", "Gene", "PTATO", "BeforePTATO","RemovedPTATO",  "Type", "CallableLoci")], 
                        id.vars = c("Label", "Days_after_sort", "Gene", "Type", "CallableLoci"))
Unique_indel_counts_m$Indels_Callable <- Unique_indel_counts_m$value / Unique_indel_counts_m$CallableLoci
Unique_indel_counts_m <- merge(Unique_indel_counts_m, IndelCounts_nonPTA_summary, by = "Gene", all.x = T)

Unique_indel_counts_m$Expected <- Unique_indel_counts_m$Indels_per_day_mean * Unique_indel_counts_m$Days_after_sort
Unique_indel_counts_m$ObsVsExp <- Unique_indel_counts_m$Indels_Callable / Unique_indel_counts_m$Expected

Unique_indel_counts_m$Indels_per_day <- Unique_indel_counts_m$Indels_Callable / Unique_indel_counts_m$Days_after_sort

Unique_indel_counts_m$variable <- as.vector(Unique_indel_counts_m$variable)
Unique_indel_counts_m$variable[Unique_indel_counts_m$variable == "BeforePTATO"] <- "Before\nPTATO"
Unique_indel_counts_m$variable[Unique_indel_counts_m$variable == "RemovedPTATO"] <- "Removed\nby PTATO"

Unique_indel_counts_m$variable <- factor(Unique_indel_counts_m$variable, levels = c("PTATO","Removed\nby PTATO", "Before\nPTATO"))


Unique_indel_counts_m$Linestart <- 1
Unique_indel_counts_PTA <- Unique_indel_counts_m[Unique_indel_counts_m$Type == "PTA",]

Unique_indel_counts_PTA$Label2 <- gsub(pattern = "-", replacement = "-\n", Unique_indel_counts_PTA$Label)
Unique_indel_counts_PTA$Label2 <- factor(Unique_indel_counts_PTA$Label2, levels = c("MSH2KO-\nPTA1", "FANCCKO-\nPTA1","WT-\nPTA2","WT-\nPTA1" ))

# Print the accuracy of mutation burden prediction
Unique_indel_counts_PTA$Accuracy <- ifelse(Unique_indel_counts_PTA$ObsVsExp < 1, Unique_indel_counts_PTA$ObsVsExp, 1/Unique_indel_counts_PTA$ObsVsExp)
mean(Unique_indel_counts_PTA$Accuracy[Unique_indel_counts_PTA$variable == "PTATO"])


ggplot(Unique_indel_counts_PTA, 
       aes(x = Label2, y = ObsVsExp, col = variable )) + 
  geom_hline(yintercept = 1, linewidth = 0.2, linetype = 2) + coord_flip() +
  
  geom_linerange(aes(ymin = Linestart, ymax = ObsVsExp), linetype = 3, linewidth = 0.2) +  
  geom_point(aes(shape  = variable, fill = variable)) +
  # scale_y_continuous(breaks = c(seq(0,30, 10), 1)) +
  # scale_shape_manual(values = c(23, 21)) +
  # scale_fill_manual(values = c( "#009E73","#D55E00")) +
  # scale_color_manual(values = c("#009E73","#D55E00")) +
   theme_classic(base_size = 6) +
  labs(fill = "Filter", color = "Filter", shape = "Filter", y = "Unique Indels\nObserved/Expected", x = "Sample") +
  theme(legend.margin=margin( l = 0, r = 0, unit='cm'),
        panel.border = element_rect(colour = "black", fill=NA, linewidth= 0.5))

ggsave(paste(Output_dir_Figure2, "Fig_2C_PMCAHH1_IndelObsVsExp.pdf", sep = ""), width = 5.5, height = 3, units = "cm")


Indel_burden <- Unique_indel_counts_PTA %>% group_by(variable) %>%
  summarise(across(c(ObsVsExp, Accuracy), list(mean = mean, sd = sd)))

ggplot(Indel_burden, aes( x= variable, y = ObsVsExp_mean)) + geom_bar(stat = "identity", width = 0.8) + 
  
  geom_errorbar(aes(ymin=ObsVsExp_mean-ObsVsExp_sd, ymax=ObsVsExp_mean+ObsVsExp_sd), width=.2,
                position=position_dodge(.9))  +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.4) +
 #geom_text(aes(y = 1.9, label = round(Accuracy_mean,2)), size = 1.8) +
  coord_flip() + 
  scale_y_continuous(expand = c(0,0.01), limits = c(0,28)) +
  theme_classic(base_size = 6) +
  labs(x = "Unique indels", y= "Indel burden\n(Observed vs Expected)")
ggsave(paste(Output_dir_Figure2, "Fig_2E_PMCAHH1_IndelObsVsExp2.pdf", sep = ""), width = 4, height = 3, units = "cm")


# Plot the indel rate per day in culture
Unique_counts_m2 <- Unique_counts_m

Unique_counts_m2 <- Unique_counts_m2[Unique_counts_m2$variable != "SCAN2" & Unique_counts_m2$variable != "Removed\nby PTATO",]

Unique_counts_m2$X_label <- gsub(Unique_counts_m2$Label, pattern = ".*-", replacement = "")
Unique_counts_m2$X_label[Unique_counts_m2$variable == "Before\nPTATO" & Unique_counts_m2$Type == "PTA"] <- paste(Unique_counts_m2$X_label[Unique_counts_m2$variable == "Before\nPTATO" & Unique_counts_m2$Type == "PTA"], " Before", sep = "")
Unique_counts_m2$X_label[Unique_counts_m2$variable != "Before\nPTATO" & Unique_counts_m2$Type == "PTA"] <- paste(Unique_counts_m2$X_label[Unique_counts_m2$variable != "Before\nPTATO" & Unique_counts_m2$Type == "PTA"], " After", sep = "")

Unique_counts_m2 <- Unique_counts_m2[-which(Unique_counts_m2$Type == "Subclone" & Unique_counts_m2$variable == "PTATO"),]
Unique_counts_m2$X_label <- factor(Unique_counts_m2$X_label, levels = c("Subclone1", "Subclone2", "PTA1 Before", "PTA1 After","PTA2 Before" ,"PTA2 After"))
Unique_counts_m2$Gene <- factor(Unique_counts_m2$Gene, levels =  c("WT", "FANCC", "MSH2"))

# Calculate the mean indel rate for the subclones and PTA samples for each gene
Unique_counts_m2 %>% group_by(Gene, Type) %>%
  summarise(across(Muts_per_day, list(mean = mean, min = min, max = max)))

basesubstitution_rate_plot <- ggplot(Unique_counts_m2, aes(x = X_label, y = Muts_per_day, fill = X_label)) + 
  geom_bar(stat = "identity") +
  facet_grid(~Gene, scales = "free", space = "free", switch = "x") + 
  geom_text(aes(label = round(Muts_per_day,0)), vjust=-0.25, size  = 1.6) +
  theme_classic(base_size = 6) +
  scale_y_continuous(expand = c(0,0), limits = c(0,max(Unique_counts_m2$Muts_per_day)*1.1)) +
  scale_fill_manual(values = c("#33a02c", "#33a02c", "#a6cee3", "#1f78b4","#a6cee3","#1f78b4")) +
  theme(strip.placement = "outside", 
        axis.text.x = element_text(angle = 45,  hjust = 1), 
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        legend.position = "none",) +
  labs(y = "Base substitution rate\n(per day in culture)", x = "Sample")

pdf(file = paste(Output_dir_FigureS5, "FigureS5A_SBS_Rate.pdf", sep = ""), width = 3, height = 2, onefile = F)
plot_adjusted_strips(basesubstitution_rate_plot, strip_colors = c(AHH1_colors[c(1,2,3)]))
dev.off()



# Plot the indel rate per day in culture
Unique_indel_counts_m2 <- Unique_indel_counts_m
Unique_indel_counts_m2 <- Unique_indel_counts_m2[Unique_indel_counts_m2$variable != "SCAN2" & Unique_indel_counts_m2$variable != "Removed\nby PTATO",]

Unique_indel_counts_m2$X_label <- gsub(Unique_indel_counts_m2$Label, pattern = ".*-", replacement = "")
Unique_indel_counts_m2$X_label[Unique_indel_counts_m2$variable == "Before\nPTATO" & Unique_indel_counts_m2$Type == "PTA"] <- paste(Unique_indel_counts_m2$X_label[Unique_indel_counts_m2$variable == "Before\nPTATO" & Unique_indel_counts_m2$Type == "PTA"], " Before", sep = "")
Unique_indel_counts_m2$X_label[Unique_indel_counts_m2$variable != "Before\nPTATO" & Unique_indel_counts_m2$Type == "PTA"] <- paste(Unique_indel_counts_m2$X_label[Unique_indel_counts_m2$variable != "Before\nPTATO" & Unique_indel_counts_m2$Type == "PTA"], " After", sep = "")

Unique_indel_counts_m2 <- Unique_indel_counts_m2[-which(Unique_indel_counts_m2$Type == "Subclone" & Unique_indel_counts_m2$variable == "PTATO"),]
Unique_indel_counts_m2$X_label <- factor(Unique_indel_counts_m2$X_label, levels = c("Subclone1", "Subclone2", "PTA1 Before", "PTA1 After","PTA2 Before" ,"PTA2 After"))

# Calculate the mean indel rate for the subclones and PTA samples for each gene
Unique_indel_counts_m2 %>% group_by(Gene, Type) %>%
  summarise(across(Indels_per_day, list(mean = mean, min = min, max = max)))
Unique_indel_counts_m2$Genelabel <- gsub("-.*", replacement = "", Unique_indel_counts_m2$Label)
Unique_indel_counts_m2$Genelabel <- factor(Unique_indel_counts_m2$Genelabel, levels =  c("WT", "FANCCKO", "MSH2KO"))

indel_rate_plot <- ggplot(Unique_indel_counts_m2, aes(x = X_label, y = Indels_per_day, fill = X_label)) + 
  geom_bar(stat = "identity") +
  facet_grid(~Genelabel, scales = "free", space = "free", switch = "x") + 
  geom_text(aes(label = round(Indels_per_day,0)), vjust=-0.25, size  = 1.6) +
  theme_classic(base_size = 6) +
  scale_y_continuous(expand = c(0,0), limits = c(0,max(Unique_indel_counts_m2$Indels_per_day)*1.1)) +
  scale_fill_manual(values = c("#33a02c", "#33a02c", "#a6cee3", "#1f78b4","#a6cee3","#1f78b4")) +
  theme(strip.placement = "outside", 
        axis.text.x = element_text(angle = 45,  hjust = 1), 
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        legend.position = "none",) +
  labs(y = "Indel rate (per day in culture)", x = "Sample")

pdf(file = paste(Output_dir_FigureS6, "FigureS6A_IndelRate.pdf", sep = ""), width = 3, height = 2, onefile = F)
plot_adjusted_strips(indel_rate_plot, strip_colors = c(AHH1_colors[c(1,2,3)]))
dev.off()


#### Plot the raw and extrapolated number of mutation per sample (for PTATO and SCAN2)
# Only select the mutations that are not shared with the (sub)clones
sSNVs_SCAN2_unique_rescue <-  lapply(sSNVs_SCAN2_unique, function(x) x[which(x$FILTER =="RESCUE"),]) 
sSNVs_SCAN2_unique_PASS <-  lapply(sSNVs_SCAN2_unique, function(x) x[which(x$FILTER =="PASS"),]) 

# SCAN2 reports the estimated mutation burden in separate files
SCAN2_burdens <- read.delim(paste(Input_dir, "PMCAHH1_SCAN2_SNV_Burden.txt", sep = ""))

SCAN2_unique_counts <- data.frame(Sample = names(sSNVs_SCAN2_unique), 
                                  PASS = lengths(sSNVs_SCAN2_unique_PASS), 
                                  FAIL = 0,
                                  RESCUE = lengths(sSNVs_SCAN2_unique_rescue) )
SCAN2_unique_counts <- merge(SCAN2_unique_counts, Metadata, by = "Sample")
SCAN2_unique_counts$Caller <- "SCAN2"
SCAN2_unique_counts <- merge(SCAN2_unique_counts, SCAN2_burdens[,c("Sample", "Extrapolated")])
sSNVs_unique_PASS <-   lapply(sSNVs_unique_raw, function(x) x[which(x$FILTER =="PASS"),]) 
sSNVs_unique_FAIL <-   lapply(sSNVs_unique_raw, function(x) x[which(x$FILTER == "FAIL"),]) 

Unique_counts <- data.frame(Sample = names(sSNVs_unique), 
                            PASS = lengths(sSNVs_unique_PASS), 
                            FAIL = lengths(sSNVs_unique_FAIL),
                            RESCUE = 0)
Unique_counts <- merge(Unique_counts, Metadata, by = "Sample")
Unique_counts <- Unique_counts[Unique_counts$Type == "PTA",]
Unique_counts$Caller <- "PTATO"
Unique_counts$Extrapolated <- Unique_counts$PASS / Unique_counts$CallableLoci

Unique_counts <- rbind(Unique_counts, SCAN2_unique_counts)

Unique_counts <- merge(Unique_counts, VariantCounts_nonPTA_summary, by = "Gene", all.x = T)
Unique_counts$Expected <- Unique_counts$Days_after_sort * Unique_counts$Muts_per_day_mean

Unique_counts_m <- melt(Unique_counts[,c("Label", "Sample", "PASS", "FAIL", "RESCUE", "Extrapolated", "Expected", "Caller")])
Unique_counts_m$TypeCounts <-  ifelse(Unique_counts_m$variable == "Extrapolated", "Extrapolated", "Uncorrected")
Unique_counts_m$TypeCounts <- factor(Unique_counts_m$TypeCounts, levels= c( "Uncorrected","Extrapolated"))

# Unique_counts_m$fill[Unique_counts_m$Caller == "SCAN2" & Unique_counts_m$variable == "PASS"] <- "SCAN2\nRaw"
# Unique_counts_m$fill[Unique_counts_m$Caller == "SCAN2" & Unique_counts_m$variable == "RESCUE"] <- "SCAN2\nRaw"
# Unique_counts_m$fill[Unique_counts_m$Caller == "SCAN2" & Unique_counts_m$variable == "Extrapolated"] <- "SCAN2\nExtrapolated"

Unique_counts_m$variable[Unique_counts_m$variable == "Extrapolated"] <- "PASS"
Unique_counts_m$variable <- factor(Unique_counts_m$variable, levels = c("FAIL", "RESCUE", "PASS", "Expected"))

#Unique_counts_m$fill <- factor(Unique_counts_m$fill, levels = c("PTATO\nRaw", "PTATO\nExtrapolated", "SCAN2\nRaw", "SCAN2\nExtrapolated"))

Unique_counts_m$Label <- factor(Unique_counts_m$Label, levels = c("WT-PTA1", "WT-PTA2", "FANCCKO-PTA1","MSH2KO-PTA1"))

Burden_plot <- ggplot(Unique_counts_m[Unique_counts_m$variable != "Expected",], aes(x = TypeCounts, y = value, fill = variable)) +
  geom_bar(stat = "identity", col = "white", linewidth = 0.2, width = 0.9) +
  ggh4x::facet_nested(.~Label + Caller, space = "free", scales = "free", switch = "x") + 
  geom_hline(data = Unique_counts_m[Unique_counts_m$variable == "Expected",], aes(yintercept = value)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,max(Unique_counts_m$value + 100))) +
  scale_fill_manual(values = c("#D55E00", "#E6AB02","#009E73")) +
  labs(x = "Sample", y = "Unique base substitutions (#)", fill = "Filter") +
  theme_classic(base_size = 6) +
  theme(strip.placement = "outside",
        strip.background = element_rect(color = "white"),
        axis.text.x = element_text(angle =45, hjust = 1),
        legend.key.size = unit(0.25, 'cm')) 
Burden_plot
pdf(file = paste(Output_dir_FigureS5, "FigureS5B_SBSburden.pdf", sep = ""), width = 4, height = 2, onefile = F)
plot_adjusted_strips(Burden_plot, strip_colors = c(AHH1_colors[c(1,1,2,3)]))
dev.off()

#### Plot the raw and extrapolated number of indels per sample (for PTATO)
# Only select the indels that are not shared with the (sub)clones

Unique_indel_counts <- data.frame(Sample = names(sIndels_unique), 
                            PASS = lengths(sIndels_unique), 
                            RECURRENT = lengths(lapply(sIndels_unique_raw, function(x) x[which(x$FILTER =="RECURRENT"),]) ),
                            HOMOPOLYMER = lengths(lapply(sIndels_unique_raw, function(x) x[which(x$FILTER =="HOMOPOLYMER"),]) ),
                            RECURRENT_HOMOPOLYMER = lengths(lapply(sIndels_unique_raw, function(x) x[which(x$FILTER =="RECURRENT;HOMOPOLYMER"),]) ))

Unique_indel_counts <- merge(Unique_indel_counts, Metadata, by = "Sample")
Unique_indel_counts <- Unique_indel_counts[Unique_indel_counts$Type == "PTA",]
Unique_indel_counts$Extrapolated <- Unique_indel_counts$PASS / Unique_indel_counts$CallableLoci

Unique_indel_counts <- merge(Unique_indel_counts, IndelCounts_nonPTA_summary, by = "Gene", all.x = T)
Unique_indel_counts$Expected <- Unique_indel_counts$Days_after_sort * Unique_indel_counts$Indels_per_day_mean

Unique_indel_counts_m <- melt(Unique_indel_counts[,c("Label", "Sample", "PASS", "RECURRENT", "HOMOPOLYMER", "RECURRENT_HOMOPOLYMER", "Extrapolated", "Expected")])
Unique_indel_counts_m$TypeCounts <-  ifelse(Unique_indel_counts_m$variable == "Extrapolated", "Extrapolated", "Uncorrected")
Unique_indel_counts_m$TypeCounts <- factor(Unique_indel_counts_m$TypeCounts, levels= c( "Uncorrected","Extrapolated"))

Unique_indel_counts_m$variable[Unique_indel_counts_m$variable == "Extrapolated"] <- "PASS"
Unique_indel_counts_m$variable <- factor(Unique_indel_counts_m$variable, levels = c("RECURRENT", "RECURRENT_HOMOPOLYMER", "HOMOPOLYMER", "PASS", "Expected"))

Unique_indel_counts_m$Label <- factor(Unique_indel_counts_m$Label, levels = c("WT-PTA1", "WT-PTA2", "FANCCKO-PTA1","MSH2KO-PTA1"))


IndelBurden_plot <- ggplot(Unique_indel_counts_m[Unique_indel_counts_m$variable != "Expected",], aes(x = TypeCounts, y = value, fill = variable)) +
  geom_bar(stat = "identity", col = "white", linewidth = 0.2, width = 0.9) +
  ggh4x::facet_nested(.~Label, space = "free", scales = "free", switch = "x") + 
  geom_hline(data = Unique_indel_counts_m[Unique_indel_counts_m$variable == "Expected",], aes(yintercept = value)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,max(Unique_indel_counts_m$value + 100))) +
  scale_fill_manual(values = c("#D55E00","#a6761d", "#E6AB02","#009E73")) +
  labs(x = "Sample", y = "Unique indels (#)", fill = "Filter") +
  theme_classic(base_size = 6) +
  theme(strip.placement = "outside",
        strip.background = element_rect(color = "white"),
        axis.text.x = element_text(angle =45, hjust = 1),
        legend.key.size = unit(0.25, 'cm')) 
IndelBurden_plot
pdf(file = paste(Output_dir_FigureS6, "FigureS6B_IndelFiltering.pdf", sep = ""), width = 4, height = 2, onefile = F)
plot_adjusted_strips(IndelBurden_plot, strip_colors = c(AHH1_colors[c(1,1,2,3)]))
dev.off()



### 96-profiles and signatures

mut_mat_raw <- mut_matrix(vcf_list = sSNVs_unique_RAW, ref_genome = ref_genome)
colnames(mut_mat_raw) <- paste(colnames(mut_mat_raw), "RAW", sep = "_")

mut_mat_unique <- mut_matrix(vcf_list = sSNVs_unique, ref_genome = ref_genome)
colnames(mut_mat_unique) <- paste(colnames(mut_mat_unique), "PTATO", sep = "_")

sSNVs_unique_fail <- lapply(sSNVs_unique_raw, function(x) x[which(x$FILTER == "FAIL"),])
mut_mat_fail <- mut_matrix(vcf_list = sSNVs_unique_fail, ref_genome = ref_genome)
colnames(mut_mat_fail) <- paste(colnames(mut_mat_fail), "FAIL", sep = "_")

mut_mat_SCAN2 <- mut_matrix(vcf_list = sSNVs_SCAN2_unique, ref_genome = ref_genome)
colnames(mut_mat_SCAN2) <-  paste(colnames(mut_mat_SCAN2), "SCAN2", sep = "_")
mut_mat_merged <- cbind(mut_mat_raw, mut_mat_unique, mut_mat_fail, mut_mat_SCAN2)

# Calculate the relative contribution of each mutation channel
mut_mat_merged2 <- sweep(mut_mat_merged,2,colSums(mut_mat_merged),`/`)
mut_mat_merged3 <- mut_mat_merged2

# Calculate the mean relative contribution of the WT PTA samples
mut_mat_merged3 <- cbind(mut_mat_merged3, rowMeans(mut_mat_merged2[,c("PMCAHH1-WT-C6SC1_PTATO", "PMCAHH1-WT-C19SC1_PTATO")]))
colnames(mut_mat_merged3)[ncol(mut_mat_merged3)] <- "WT_PTATO"

mut_mat_merged3 <- cbind(mut_mat_merged3, rowMeans(mut_mat_merged2[,c("PMCAHH1-WT-C6SC1_RAW", "PMCAHH1-WT-C19SC1_RAW")]))
colnames(mut_mat_merged3)[ncol(mut_mat_merged3)] <- "WT_RAW"

mut_mat_merged3 <- cbind(mut_mat_merged3, rowMeans(mut_mat_merged2[,c("PMCAHH1-WT-C6SC1_FAIL", "PMCAHH1-WT-C19SC1_FAIL")]))
colnames(mut_mat_merged3)[ncol(mut_mat_merged3)] <- "WT_FAIL"

mut_mat_merged3 <- cbind(mut_mat_merged3, rowMeans(mut_mat_merged2[,c("PMCAHH1-WT-C6SC1_SCAN2", "PMCAHH1-WT-C19SC1_SCAN2")]))
colnames(mut_mat_merged3)[ncol(mut_mat_merged3)] <- "WT_SCAN2"

# Signatures
signatures <- get_known_signatures()
signatures <- signatures[,1:30]
pta_signatures <- read.table(paste(Input_dir, "Resources/PTA_Artefact_Signature.txt", sep = ""), header=T)
signatures <- cbind(signatures,pta_signatures$PTA)
colnames(signatures)[length(colnames(signatures))] <- "PTA"
merged_signatures <- merge_signatures(signatures, cos_sim_cutoff = 0.85)


mut_mat2 <- cbind(mut_mat_merged3, pta_signatures)
cos_sim_samples_signatures <- cos_sim_matrix(mut_mat2, mut_mat2)

cos_sim_samples_signatures_m <- melt(cos_sim_samples_signatures)

cos_sim_samples_signatures_m$Type <- gsub(pattern = ".*_", replacement = "", cos_sim_samples_signatures_m$Var1)
cos_sim_samples_signatures_m$X <-  gsub(pattern = "_.*", replacement = "", cos_sim_samples_signatures_m$Var2)
cos_sim_samples_signatures_m$ID <-  gsub(pattern = "_.*", replacement = "", cos_sim_samples_signatures_m$Var1)

cos_sim_samples_signatures_m <- cos_sim_samples_signatures_m[cos_sim_samples_signatures_m$Var2 %in% c("PTA", "PMCAHH1-MSH2KO-C27E06SC51B06_RAW", "PMCAHH1-FANCCKO-C02B03SC03E05_RAW", "AHH1WTG2SCG6_RAW"),]

cos_sim_samples_signatures_m <- cos_sim_samples_signatures_m[cos_sim_samples_signatures_m$ID %in% c("WT",
                                                                                                    "PMCAHH1-FANCCKO-C02B03SC03E05-PTAP1D7",
                                                                                                    "PMCAHH1-MSH2KO-C27E06SC51B06-PTAP1E7"),]

cos_sim_samples_signatures_m[cos_sim_samples_signatures_m$ID == "PMCAHH1-FANCCKO-C02B03SC03E05-PTAP1D7" & cos_sim_samples_signatures_m$X == "PMCAHH1-MSH2KO-C27E06SC51B06",]

# Adjust the sample labels for plotting
cos_sim_samples_signatures_m$Sample <- cos_sim_samples_signatures_m$ID
cos_sim_samples_signatures_m$Sample[cos_sim_samples_signatures_m$ID ==  "PMCAHH1-FANCCKO-C02B03SC03E05-PTAP1D7"] <- "FANCCKO-PTA"
cos_sim_samples_signatures_m$Sample[cos_sim_samples_signatures_m$ID ==  "PMCAHH1-MSH2KO-C27E06SC51B06-PTAP1E7"] <- "MSH2KO-PTA"
cos_sim_samples_signatures_m$Sample[cos_sim_samples_signatures_m$ID ==  "WT"] <- "WT-PTA"
cos_sim_samples_signatures_m$Sample <- factor(cos_sim_samples_signatures_m$Sample, levels = c("WT-PTA","FANCCKO-PTA", "MSH2KO-PTA" ))

cos_sim_samples_signatures_m$X_labels <- cos_sim_samples_signatures_m$X
cos_sim_samples_signatures_m$X_labels[cos_sim_samples_signatures_m$X == "AHH1WTG2SCG6"] <- "WT-\nSubclone1"
cos_sim_samples_signatures_m$X_labels[cos_sim_samples_signatures_m$X == "PMCAHH1-FANCCKO-C02B03SC03E05"] <- "FANCCKO-\nSubclone1"
cos_sim_samples_signatures_m$X_labels[cos_sim_samples_signatures_m$X == "PMCAHH1-MSH2KO-C27E06SC51B06"] <- "MSH2KO-\nSubclone2"
cos_sim_samples_signatures_m$X_labels <- factor(cos_sim_samples_signatures_m$X_labels, levels = c("WT-\nSubclone1",  "FANCCKO-\nSubclone1","MSH2KO-\nSubclone2", "PTA"))
# Reorder the Y-axis labels
cos_sim_samples_signatures_m$Type <- factor(cos_sim_samples_signatures_m$Type, levels = c("SCAN2", "PTATO", "FAIL", "RAW"))

# plot_labels <- data.frame(ID = c(Inputs$ID, "WT", "PTA", "AHH1WTG2SCG6"), 
#                           X = c(Inputs$ID, "WT", "PTA", "AHH1WTG2SCG6"), 
#                           Label_X = c("MSH2KO-\nSC1", "MSH2KO-\nSC2", "MSH2KO-\nPTA", "FANCCKO-\nSC1",
#                                       "FANCCKO-\nPTA", "WT-\nPTA1", "WT-\nPTA2", "WT-\nPTA3","WT-\nPTA", "PTA", "WT-\nSC1"),
#                           Label_Y = c("MSH2KO-SC1", "MSH2KO-SC2", "MSH2KO-PTA", "FANCCKO-SC1",
#                                       "FANCCKO-PTA", "WT-PTA1", "WT-PTA2", "WT-PTA3","WT-PTA", "PTA","WT-SC1"))
# cos_sim_samples_signatures_m <- merge(cos_sim_samples_signatures_m, plot_labels[,c("ID", "Label_Y")], by = "ID", all.x = T)
# cos_sim_samples_signatures_m <- merge(cos_sim_samples_signatures_m, plot_labels[,c("X", "Label_X")], by = "X", all.x = T)
# cos_sim_samples_signatures_m$Label_X <- factor(cos_sim_samples_signatures_m$Label_X, levels = c("WT-\nSC1", "FANCCKO-\nSC1", "MSH2KO-\nSC2","PTA" ))

fig2c_raw <- ggplot(cos_sim_samples_signatures_m, aes(x = X_labels, y = Type, fill = value, label = round(value,2))) + geom_raster() +
  geom_text(size = 1.8) + 
  facet_grid(Sample~., switch = "both", scales = "free") +
  scale_fill_gradient(low = "white", high = "firebrick3") +
  scale_x_discrete(position = "top", expand = c(0,0)) +
  theme_classic(base_size = 6) +
  theme( strip.placement = "outside", 
         axis.text.x = element_text(angle = -45, hjust = 1),
         legend.position = "bottom") +
  labs(x = "Base Substitution Profile", y= "Sample", fill = "Cosine\nSimilarity") + 
  guides(fill = guide_colourbar(barwidth = 5, barheight = 0.5))

fig2c_raw

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

# dev.off()
# pdf(file = paste(Output_dir_Figure2, "Fig_2d_PMCAHH1_CosineHeatmap_SBS.pdf", sep = ""), width = 1.5, height = 2.8, onefile = F)
# plot_adjusted_strips(fig2c_raw, strip_colors = AHH1_colors[c(2,3, 1)])
# dev.off()


### 
mut_mat2 <- cbind(mut_mat_merged2, pta_signatures)
cos_sim_samples_signatures <- cos_sim_matrix(mut_mat2, mut_mat2)

cos_sim_samples_signatures_m <- melt(cos_sim_samples_signatures)

cos_sim_samples_signatures_m$Type <- gsub(pattern = ".*_", replacement = "", cos_sim_samples_signatures_m$Var1)
cos_sim_samples_signatures_m$X <-  gsub(pattern = "_.*", replacement = "", cos_sim_samples_signatures_m$Var2)
cos_sim_samples_signatures_m$ID <-  gsub(pattern = "_.*", replacement = "", cos_sim_samples_signatures_m$Var1)

cos_sim_samples_signatures_m <- cos_sim_samples_signatures_m[cos_sim_samples_signatures_m$Var2 %in% c("PTA", "PMCAHH1-MSH2KO-C27E06SC51B06_RAW", "PMCAHH1-FANCCKO-C02B03SC03E05_RAW", "AHH1WTG2SCG6_RAW"),]

cos_sim_samples_signatures_m <- cos_sim_samples_signatures_m[cos_sim_samples_signatures_m$ID %in% c("PMCAHH1-WT-C19SC1", "PMCAHH1-WT-C6SC1",
                                                                                                    "PMCAHH1-FANCCKO-C02B03SC03E05-PTAP1D7",
                                                                                                    "PMCAHH1-MSH2KO-C27E06SC51B06-PTAP1E7"),]

# Keep only the cosine similarities between the PTA sample and the corresponding bulk sample (so remove the cos sims between PTA samples and the spectra of the other cell lines)
cos_sim_overview <- rbind(cos_sim_samples_signatures_m[cos_sim_samples_signatures_m$Var2 == "AHH1WTG2SCG6_RAW" & cos_sim_samples_signatures_m$ID %in% c("PMCAHH1-WT-C19SC1", "PMCAHH1-WT-C6SC1"),],
                                      cos_sim_samples_signatures_m[cos_sim_samples_signatures_m$Var2 == "PMCAHH1-FANCCKO-C02B03SC03E05_RAW" & cos_sim_samples_signatures_m$ID %in% c("PMCAHH1-FANCCKO-C02B03SC03E05-PTAP1D7"),],
                                      cos_sim_samples_signatures_m[cos_sim_samples_signatures_m$Var2 == "PMCAHH1-MSH2KO-C27E06SC51B06_RAW" & cos_sim_samples_signatures_m$ID %in% c("PMCAHH1-MSH2KO-C27E06SC51B06-PTAP1E7"),],
                                      cos_sim_samples_signatures_m[cos_sim_samples_signatures_m$Var2 == "PTA",])
cos_sim_overview$Var2 <- as.vector(cos_sim_overview$Var2)
cos_sim_overview$Var2[cos_sim_overview$Var2 != "PTA"] <- "Subclone"

cos_sim_means <- cos_sim_overview %>% group_by(Var2, Type) %>% summarise(across(value, list(mean = mean, sd = sd)))
cos_sim_means$Var2 <- factor(cos_sim_means$Var2, levels = c("Subclone", "PTA"))
cos_sim_means$Label <- cos_sim_means$Type
cos_sim_means$Label[cos_sim_means$Label == "FAIL"] <- "RemovedPTATO"
cos_sim_means$Label[cos_sim_means$Label == "RAW"] <- "BeforePTATO"
cos_sim_means$Label <- factor(cos_sim_means$Label, levels = c("SCAN2", "PTATO", "RemovedPTATO", "BeforePTATO"))

ggplot(cos_sim_means, aes(x = Var2, y = Label, fill = value_mean, label = round(value_mean,2))) + 
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(size = 1.8) + 
  scale_fill_gradient(low = "white", high = "firebrick3") +
  #scale_x_discrete(position = "top") +
  #scale_x_discrete(position = "top", expand = c(0,0)) +
  theme_classic(base_size = 6) +
  theme(plot.title = element_text(hjust = 0.5, size = 6),
         legend.position = "bottom" ) + ggtitle("Cosine similarity") +
  labs(x = "Mutational profile", y= "Sample", fill = "Cosine\nSimilarity") + 
  guides(fill = guide_colourbar(barwidth = 2.5, barheight = 0.5, ))

ggsave(paste(Output_dir_Figure2, "Fig_2F_PMCAHH1_Cossim.pdf", sep = ""), width = 4, height = 4.25, units = "cm")


### REFIT
mut_mat_merged2
contri_boots <- data.frame()
for(Sample in Metadata$Sample[Metadata$Type == "PTA"]){
  print(Sample)
  
  Subclone_Sample <-  Metadata$Sample[Metadata$Type == "Subclone" & Metadata$Gene == Metadata$Gene[Metadata$Sample == Sample]][length(Metadata$Sample[Metadata$Type == "Subclone" & Metadata$Gene == Metadata$Gene[Metadata$Sample == Sample]])]
  
  mut_mats_sample <- mut_mat_merged[,grep(Sample, x = colnames(mut_mat_merged))]
  
  Signatures_sample <- cbind(PTA = pta_signatures[,1], Subclone = as.numeric(mut_mat_merged2[,paste(Subclone_Sample, "RAW", sep = "_")]))
  
  contri_boots_sample <- fit_to_signatures_bootstrapped(mut_mats_sample,
                                                        Signatures_sample,
                                                 n_boots = 100,
                                                 method = "strict", max_delta = 0.005)
  
  contri_boots <- rbind(contri_boots, contri_boots_sample)
}

contri_boots_rel <- contri_boots / rowSums(contri_boots)

# Change to long format for plotting
contri_tb <- contri_boots_rel %>% as.data.frame() %>% tibble::rownames_to_column("exp") %>% 
  tidyr::gather(key = "sig", value = "contri", -exp) %>% 
  dplyr::mutate(Sample = gsub("_[^_]+$", "", exp), Sample = factor(Sample, 
                                                                   levels = unique(Sample)), sig = factor(sig, levels = unique(sig)))
# Calculate the mean contribution for each signature
contri_tb2 <- contri_tb %>% dplyr::group_by(Sample, sig) %>% 
  dplyr::summarise(mean = mean(contri)) %>% 
  dplyr::ungroup()

contri_tb2$Type <- gsub(contri_tb2$Sample, pattern = ".*_", replacement = "")

# Calculate the mean contribution for each signature
contri_tb_means <- contri_tb2 %>% dplyr::group_by(Type, sig) %>% 
  dplyr::summarise(across(mean, list(mean = mean, sd = sd)))
contri_tb_means$Label <- contri_tb_means$Type
contri_tb_means$Label[contri_tb_means$Label == "FAIL"] <- "RemovedPTATO"
contri_tb_means$Label[contri_tb_means$Label == "RAW"] <- "BeforePTATO"
contri_tb_means$Label <- factor(contri_tb_means$Label, levels = c("SCAN2", "PTATO", "RemovedPTATO", "BeforePTATO"))
ggplot(contri_tb_means, aes(x = Label, y  = mean_mean, fill = sig)) + geom_bar(stat = "identity", width = 0.8) +   
  geom_errorbar(data = contri_tb_means[contri_tb_means$sig == "Subclone",], aes(ymin=mean_mean-mean_sd, ymax=mean_mean+mean_sd), width=.2)  +
  geom_text(data = contri_tb_means[contri_tb_means$sig == "Subclone",], aes(label = round(mean_mean,2), y = 1.1), size = 1.8, fontface = 3) +
  coord_flip() + 
  scale_y_continuous(breaks = seq(0,1, 0.25), labels = scales::percent)+
  scale_fill_manual(values = c( "#D55E00", "#009E73")) +
  theme_classic(base_size = 6) +
  labs(y = "Base substitutions (%)", fill = "Mutational profile", title = "Strict signature refit") +
  theme(legend.key.size = unit(0.25, 'cm'), legend.position = "bottom", plot.title = element_text(hjust = 0.5, size = 6)) +
  guides(fill = guide_legend(title.position = "top", 
                              title.hjust = 0.5) )

ggsave(paste(Output_dir_Figure2, "Fig_2H_PMCAHH1_Refit.pdf", sep = ""), width = 5, height = 4.3, units = "cm")


  ### Indels

sIndels_unique_raw <- get_mut_type(sIndels_unique_raw, type = "indel")
sIndels_unique_raw <- get_indel_context(sIndels_unique_raw, ref_genome = ref_genome)
indel_counts_raw <- count_indel_contexts(sIndels_unique_raw)

indel_counts_fail <- count_indel_contexts( lapply(sIndels_unique_raw, function(x) x[which(x$FILTER != "PASS"),]) )

## Takes a long time
# exclusion_context <- GRanges()
# for(chr in seqlevels(indel_recurrency_grl)){
#   print(chr)
#   exclusion_chr <- indel_recurrency_grl[which(seqnames(indel_recurrency_grl) == chr)]
#   exclusion_chr_context <- get_indel_context(exclusion_chr, ref_genome)
#   exclusion_context <- c(exclusion_context, exclusion_chr_context)
# }
# 
# saveRDS(exclusion_context, file = paste(Input_dir, "Indels/IndelExclusionList.rds", sep = ""))

exclusion_context <- readRDS(file = paste(Input_dir, "Indels/IndelExclusionList.rds", sep = ""))

# exclusion_counts <- count_indel_contexts(exclusion_context)
# colnames(exclusion_counts) <- "PTATO\nIndel\nExclusion\nList"

indel_counts_recurrent <- count_indel_contexts(exclusion_context)
colnames(indel_counts_recurrent) <- "Recurrency\nList"
sIndels_unique <- get_mut_type(sIndels_unique, type = "indel")
sIndels_unique <- get_indel_context(sIndels_unique, ref_genome = ref_genome)
indel_counts <- count_indel_contexts(sIndels_unique)

#plot_indel_contexts(indel_counts_fail, condensed = TRUE)

cos_sim_indels_raw <- cos_sim_matrix(indel_counts_raw, cbind(indel_counts[,c("PMCAHH1-FANCCKO-C02B03SC03E05", "PMCAHH1-MSH2KO-C27E06SC51B06", "AHH1WTG2SCG6")], 
                                                             indel_counts_recurrent))
cos_sim_indels_raw <- rbind(cos_sim_indels_raw, colMeans(cos_sim_indels_raw[row.names(cos_sim_indels_raw) %in% c("PMCAHH1-WT-C19SC1", "PMCAHH1-WT-C6SC1"),]))
rownames(cos_sim_indels_raw)[nrow(cos_sim_indels_raw)] <- "WT"

cos_sim_indels_fail <- cos_sim_matrix(indel_counts_fail, cbind(indel_counts[,c("PMCAHH1-FANCCKO-C02B03SC03E05", "PMCAHH1-MSH2KO-C27E06SC51B06", "AHH1WTG2SCG6")], 
                                                             indel_counts_recurrent))
cos_sim_indels_fail <- rbind(cos_sim_indels_fail, colMeans(cos_sim_indels_fail[row.names(cos_sim_indels_fail) %in% c("PMCAHH1-WT-C19SC1", "PMCAHH1-WT-C6SC1"),]))
rownames(cos_sim_indels_fail)[nrow(cos_sim_indels_fail)] <- "WT"

cos_sim_indels <- cos_sim_matrix(indel_counts, cbind(indel_counts[,c("PMCAHH1-FANCCKO-C02B03SC03E05", "PMCAHH1-MSH2KO-C27E06SC51B06", "AHH1WTG2SCG6")], 
                                 indel_counts_recurrent))
cos_sim_indels <- rbind(cos_sim_indels, colMeans(cos_sim_indels[row.names(cos_sim_indels) %in% c("PMCAHH1-WT-C19SC1", "PMCAHH1-WT-C6SC1"),]))
rownames(cos_sim_indels)[nrow(cos_sim_indels)] <- "WT"

cos_sim_indels_raw_m <- melt(cos_sim_indels_raw)
cos_sim_indels_raw_m$Filter <- "RAW"

cos_sim_indels_fail_m <- melt(cos_sim_indels_fail)
cos_sim_indels_fail_m$Filter <- "FAIL"

cos_sim_indels_m <- melt(cos_sim_indels)
cos_sim_indels_m$Filter <- "PTATO"

cos_sim_indels_merged <- rbind(cos_sim_indels_raw_m,cos_sim_indels_fail_m, cos_sim_indels_m)
colnames(cos_sim_indels_merged) <- c("Sample", "Profile", "CosineSimilarity", "Filter")
cos_sim_indels_merged$Sample <- as.character(cos_sim_indels_merged$Sample)
cos_sim_indels_merged$Profile <- as.character(cos_sim_indels_merged$Profile)

# Adjust the sample labels for plotting
cos_sim_indels_merged$Sample[cos_sim_indels_merged$Sample ==  "PMCAHH1-FANCCKO-C02B03SC03E05-PTAP1D7"] <- "FANCCKO-PTA"
cos_sim_indels_merged$Sample[cos_sim_indels_merged$Sample ==  "PMCAHH1-MSH2KO-C27E06SC51B06-PTAP1E7"] <- "MSH2KO-PTA"
cos_sim_indels_merged$Sample[cos_sim_indels_merged$Sample ==  "WT"] <- "WT-PTA"
cos_sim_indels_merged$Sample <- factor(cos_sim_indels_merged$Sample, levels = c("WT-PTA","FANCCKO-PTA", "MSH2KO-PTA" ))

cos_sim_indels_merged$Profile[cos_sim_indels_merged$Profile == "AHH1WTG2SCG6"] <- "WT-\nSubclone1"
cos_sim_indels_merged$Profile[cos_sim_indels_merged$Profile == "PMCAHH1-FANCCKO-C02B03SC03E05"] <- "FANCCKO-\nSubclone1"
cos_sim_indels_merged$Profile[cos_sim_indels_merged$Profile  == "PMCAHH1-MSH2KO-C27E06SC51B06"] <- "MSH2KO-\nSubclone2"
cos_sim_indels_merged$Profile <- factor(cos_sim_indels_merged$Profile, levels = c("WT-\nSubclone1",  "FANCCKO-\nSubclone1","MSH2KO-\nSubclone2", "Recurrency\nList"))

cos_sim_indels_merged <- cos_sim_indels_merged[!is.na(cos_sim_indels_merged$Sample),]

# Reorder the Y-axis labels
cos_sim_indels_merged$Filter <- factor(cos_sim_indels_merged$Filter, levels = c("PTATO", "FAIL", "RAW"))

indel_cossim_plot <- ggplot(cos_sim_indels_merged, aes(x = Profile, y = Filter, fill = CosineSimilarity, label = round(CosineSimilarity,2))) +
  geom_raster() +
  geom_text(size = 1.8) + 
  facet_grid(Sample~., switch = "both", scales = "free") +
  scale_fill_gradient(low = "white", high = "dodgerblue4") +
  scale_x_discrete(position = "top", expand = c(0,0)) +
  theme_classic(base_size = 6) +
  theme( strip.placement = "outside", axis.text.x = element_text(angle = -45, hjust = 1), legend.position = "bottom") +
  labs(x = "Indel Profile", y= "Sample", fill = "Cosine\nSimilarity") + 
  guides(fill = guide_colourbar(barwidth = 5, barheight = 0.5))

# Old version
indel_cossim_plot

# dev.off()
# pdf(file = paste(Output_dir_Figure2, "Fig_2G_PMCAHH1_Indel_CosineHeatmap.pdf", sep = ""), width = 1.5, height = 2.8, onefile = F)
# plot_adjusted_strips(indel_cossim_plot, strip_colors = AHH1_colors[c(2,3, 1)])
# dev.off()

### 

cos_sim_indels_merged <- rbind(cos_sim_indels_raw_m,cos_sim_indels_fail_m, cos_sim_indels_m)
colnames(cos_sim_indels_merged) <- c("Sample", "Profile", "CosineSimilarity", "Filter")
cos_sim_indels_merged$Sample <- as.character(cos_sim_indels_merged$Sample)
cos_sim_indels_merged$Profile <- as.character(cos_sim_indels_merged$Profile)

cos_sim_indels_merged <- cos_sim_indels_merged[cos_sim_indels_merged$Sample %in% c("PMCAHH1-WT-C19SC1", "PMCAHH1-WT-C6SC1", "PMCAHH1-FANCCKO-C02B03SC03E05-PTAP1D7", "PMCAHH1-MSH2KO-C27E06SC51B06-PTAP1E7"),]
# Keep only the cosine similarities between the PTA sample and the corresponding bulk sample (so remove the cos sims between PTA samples and the spectra of the other cell lines)
cos_sim_indel_overview <- rbind(cos_sim_indels_merged[cos_sim_indels_merged$Profile == "AHH1WTG2SCG6" & cos_sim_indels_merged$Sample %in% c("PMCAHH1-WT-C19SC1", "PMCAHH1-WT-C6SC1"),],
                          cos_sim_indels_merged[cos_sim_indels_merged$Profile == "PMCAHH1-FANCCKO-C02B03SC03E05" & cos_sim_indels_merged$Sample %in% c("PMCAHH1-FANCCKO-C02B03SC03E05-PTAP1D7"),],
                          cos_sim_indels_merged[cos_sim_indels_merged$Profile == "PMCAHH1-MSH2KO-C27E06SC51B06" & cos_sim_indels_merged$Sample %in% c("PMCAHH1-MSH2KO-C27E06SC51B06-PTAP1E7"),],
                          cos_sim_indels_merged[cos_sim_indels_merged$Profile == "Recurrency\nList",])
cos_sim_indel_overview$Profile <- as.vector(cos_sim_indel_overview$Profile)
cos_sim_indel_overview$Profile[cos_sim_indel_overview$Profile != "Recurrency\nList"] <- "Subclone"

cos_sim_indels_means <- cos_sim_indel_overview %>% group_by(Profile, Filter) %>% summarise(across(CosineSimilarity, list(mean = mean, sd = sd)))
cos_sim_indels_means$Profile <- factor(cos_sim_indels_means$Profile, levels = c("Subclone", "Recurrency\nList"))
cos_sim_indels_means$Label <- cos_sim_indels_means$Filter
cos_sim_indels_means$Label[cos_sim_indels_means$Label == "FAIL"] <- "RemovedPTATO"
cos_sim_indels_means$Label[cos_sim_indels_means$Label == "RAW"] <- "BeforePTATO"
cos_sim_indels_means$Label <- factor(cos_sim_indels_means$Label, levels = c("SCAN2", "PTATO", "RemovedPTATO", "BeforePTATO"))

ggplot(cos_sim_indels_means, aes(x = Profile, y = Label, fill = CosineSimilarity_mean, label = round(CosineSimilarity_mean,2))) + 
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(size = 1.8) + 
  scale_fill_gradient(low = "white", high = "dodgerblue4") +
  #scale_x_discrete(position = "top") +
  #scale_x_discrete(position = "top", expand = c(0,0)) +
  theme_classic(base_size = 6) +
  theme(plot.title = element_text(hjust = 0.5, size = 6) ) + ggtitle("Cosine similarity") +
  labs(x = "Mutational profile", y= "Sample", fill = "Cosine\nSimilarity") + 
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 2.5))

ggsave(paste(Output_dir_Figure2, "Fig_2G_PMCAHH1_Indel_Cossim.pdf", sep = ""), width = 5.25, height = 3.25, units = "cm")



### 96-profiles
plot_96_profile_small <- function(mut_matrix, colors = NA, ymax = 0.2, condensed = FALSE) {
  
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
  
  # Make contribution relative
  norm_mut_matrix <- apply(mut_matrix, 2, function(x) x / sum(x))
  
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
  
  return(plot)
}

mut_mat_WT <- mut_mat_merged3[,c("AHH1WTG2SCG6_RAW", "WT_RAW", "WT_FAIL", "WT_PTATO", "WT_SCAN2" )]
colnames(mut_mat_WT) <- c("WT-\nSubclone1", "WT-PTA\nRAW", "WT-PTA\nFAIL", "WT-PTA\nPTATO", "WT-PTA\nSCAN2")

plotWT <- plot_96_profile_small(mut_mat_WT, condensed = TRUE, ymax = 0.1)
pdf(file = paste(Output_dir_FigureS5, "Figure_S5C_96profile_WT.pdf", sep = ""), width = 3.2, height = 2.5, onefile = F)

plot_adjusted_strips(plotWT, strip_colors = c(
  "#2EBAED", "#000000", "#DE1C14",
  "#D4D2D2", "#ADCC54", "#F0D0CE"
))
dev.off()


mut_mat_FANCC <- mut_mat_merged3[,c("PMCAHH1-FANCCKO-C02B03SC03E05_RAW",
                                    "PMCAHH1-FANCCKO-C02B03SC03E05-PTAP1D7_RAW",
                                    "PMCAHH1-FANCCKO-C02B03SC03E05-PTAP1D7_FAIL",
                                    "PMCAHH1-FANCCKO-C02B03SC03E05-PTAP1D7_PTATO", 
                                    "PMCAHH1-FANCCKO-C02B03SC03E05-PTAP1D7_SCAN2")]
colnames(mut_mat_FANCC) <- c("FANCCKO-\nSubclone1", "FANCCKO-\nPTA_RAW", "FANCCKO-\nPTA_FAIL", "FANCCKO-\nPTA_PTATO", "FANCCKO-\nPTA_SCAN2")

plotFANCC <- plot_96_profile_small(mut_mat_FANCC, condensed = TRUE, ymax = 0.1)
pdf(file = paste(Output_dir_FigureS5, "Figure_S5D_96profile_FANCC.pdf", sep = ""), width = 3.2, height = 2.5, onefile = F)

plot_adjusted_strips(plotFANCC, strip_colors = c(
  "#2EBAED", "#000000", "#DE1C14",
  "#D4D2D2", "#ADCC54", "#F0D0CE"
))
dev.off()

mut_mat_MSH2 <- mut_mat_merged3[,c("PMCAHH1-MSH2KO-C27E06SC51B06_RAW", 
                                   "PMCAHH1-MSH2KO-C27E06SC51B06-PTAP1E7_RAW", 
                                   "PMCAHH1-MSH2KO-C27E06SC51B06-PTAP1E7_FAIL", 
                                   "PMCAHH1-MSH2KO-C27E06SC51B06-PTAP1E7_PTATO", 
                                   "PMCAHH1-MSH2KO-C27E06SC51B06-PTAP1E7_SCAN2")]
colnames(mut_mat_MSH2) <- c("MSH2KO-\nSubclone1", "MSH2KO-\nPTA_RAW", "MSH2KO-\nPTA_FAIL", "MSH2KO-\nPTA_PTATO", 
                            "MSH2KO-\nPTA_SCAN2")

plotMSH2 <- plot_96_profile_small(mut_mat_MSH2, condensed = TRUE, ymax = 0.1)
pdf(file = paste(Output_dir_FigureS5, "Figure_S5E_96profile_MSH2.pdf", sep = ""), width = 3.2, height = 2.5, onefile = F)

plot_adjusted_strips(plotMSH2, strip_colors = c(
  "#2EBAED", "#000000", "#DE1C14",
  "#D4D2D2", "#ADCC54", "#F0D0CE"
))
dev.off()

### Signature refit
# Signatures
signatures <- get_known_signatures()
#signatures <- signatures[,1:36]
pta_signature <- read.table(paste(Input_dir, "Resources/PTA_Artefact_Signature.txt", sep = ""), header=T)
signatures <- cbind(signatures,pta_signature$PTA)
colnames(signatures)[length(colnames(signatures))] <- "PTA"
merged_signatures <- merge_signatures(signatures, cos_sim_cutoff = 0.8)

mut_mat_combined <- cbind(mut_mat_WT, mut_mat_FANCC, mut_mat_MSH2)

nmf_res <- extract_signatures(mut_mat_combined, rank = 3, nrun = 100, single_core = TRUE)
nmf_res <- rename_nmf_signatures(nmf_res, merged_signatures, cutoff = 0.75)
colnames(nmf_res$signatures) <- paste(c("A","B", "C"), "\n", colnames(nmf_res$signatures), sep = "")

pdf(file = paste(Output_dir_FigureS5, "Figure_S5F_NMF_Signatures.pdf", sep = ""), width = 4, height = 2, onefile = F)
plot_96_profile_small(nmf_res$signatures, condensed = TRUE, ymax = 0.1)
dev.off()

# plot_contribution(nmf_res$contribution, nmf_res$signature,
#                   mode = "relative"
# )
# 
# plot_original_vs_reconstructed(mut_mat_merged3[,!is.na(colSums(mut_mat_merged3))], nmf_res$reconstructed, 
#                                y_intercept = 0.95)


tb <- nmf_res$contribution %>% as.data.frame() %>% tibble::rownames_to_column("Signature") %>% 
  tidyr::pivot_longer(-Signature, names_to = "ID", 
                      values_to = "Contribution") %>% dplyr::mutate(ID = factor(ID, 
                                                                                    levels = unique(ID)), Signature = factor(Signature, 
                                                                                                                                 levels = unique(Signature)))
tb$Cell_line <- gsub(tb$ID, pattern = "-.*",  replacement = "")
tb$Cell_line <- factor(tb$Cell_line, levels = c("WT", "FANCCKO", "MSH2KO"))
tb$Type <- gsub(tb$ID, pattern = ".*-",  replacement = "")
tb$Type <- gsub(tb$Type, pattern = "A\n", replacement = "A_")
tb$Type <- gsub(tb$Type, pattern = "\n", replacement = "")
tb$Type <- factor(tb$Type, levels = rev(c("Subclone1", "PTA_RAW", "PTA_FAIL", "PTA_PTATO", "PTA_SCAN2")))

nmf_plot <- ggplot(tb, aes(x = Type, y = Contribution, fill = Signature)) +
  geom_bar(stat = "identity", position = position_fill(), color = "black", linewidth = 0.15, width = 0.8) +
  coord_flip()+
  facet_grid(Cell_line~., space = "free", scales = "free", switch = "y") +
  theme_classic(base_size = 6) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("#377eb8", "#e41a1c", "#4daf4a")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.placement = "outside",
        legend.key.size = unit(0.25, 'cm'))  +
  labs(x = "Sample")

pdf(file =  paste(Output_dir_FigureS5, "Figure_S5G_NMF_Contribution.pdf", sep = ""), width = 3.5, height = 2, onefile = F)
plot_adjusted_strips(nmf_plot, strip_colors = c(AHH1_colors[c(1,2,3)]))
dev.off()

sigs <- cbind(WT = mut_mat_WT[,1], PTA = pta_signature[,1])
sigs <- cbind(Subclone = mut_mat_FANCC[,1], PTA = pta_signature[,1])
sigs <- cbind(Subclone = mut_mat_MSH2[,1], PTA = pta_signature[,1])

fit_res <- fit_to_signatures(mut_mat_WT, sigs)
plot_contribution(fit_res$contribution,
                  coord_flip = FALSE,
                  mode = "absolute"
)
refit_WT <- fit_to_signatures_strict(mut_mat_WT, sigs)
refit_WT <- fit_to_signatures_strict(mut_mat_FANCC, sigs)
refit_WT <- fit_to_signatures_strict(mut_mat_MSH2, sigs)

fit_res_strict <- refit_WT$fit_res
plot_contribution(fit_res_strict$contribution,
                  coord_flip = FALSE,
                  mode = "relative"
)

### Indel profiles

indel_counts_WT <- cbind(`WT-Subclone1` = indel_counts[,"AHH1WTG2SCG6"], 
                         `WT-PTA\nRaw` = rowSums(indel_counts_raw[,c("PMCAHH1-WT-C19SC1", "PMCAHH1-WT-C6SC1")]), 
                         `WT-PTA\nPTATO` = rowSums(indel_counts[,c("PMCAHH1-WT-C19SC1", "PMCAHH1-WT-C6SC1")]))

plot_indel_contexts2(counts =indel_counts_WT, condensed = T, y_label = " (#)")

indel_counts_FANCC <- cbind(`FANCCKO-\nSubclone1` = indel_counts[,"PMCAHH1-FANCCKO-C02B03SC03E05"], 
                         `FANCCKO-PTA\nRaw` = indel_counts_raw[,"PMCAHH1-FANCCKO-C02B03SC03E05-PTAP1D7"], 
                         `FANCCKO-PTA\nPTATO` = indel_counts[,"PMCAHH1-FANCCKO-C02B03SC03E05-PTAP1D7"])

plot_indel_contexts2(counts =indel_counts_FANCC, condensed = T, y_label = " (#)") 

indel_counts_MSH2 <- cbind(`MSH2KO-\nSubclone1` = indel_counts[,"PMCAHH1-MSH2KO-C27E06SC51B06"], 
                            `MSH2KO-PTA\nRaw` = indel_counts_raw[,"PMCAHH1-MSH2KO-C27E06SC51B06-PTAP1E7"], 
                            `MSH2KO-PTA\nPTATO` = indel_counts[,"PMCAHH1-MSH2KO-C27E06SC51B06-PTAP1E7"])

plot_indel_contexts2(counts =indel_counts_MSH2,condensed = T, y_label = " (#)")

indel_counts_merged <- cbind(indel_counts_WT, indel_counts_FANCC, indel_counts_MSH2)

plot_main_indel_contexts_stacked(indel_counts_merged)

indel_counts_long <- indel_counts_merged %>% as.data.frame() %>% tibble::rownames_to_column("muttype_total") %>% 
  tidyr::separate(muttype_total, c("muttype", "muttype_sub"), 
                  sep = "_(?=[0-9])") %>% dplyr::mutate(muttype = factor(muttype, 
                                                                         levels = unique(muttype)))
indel_counts_main <- indel_counts_long %>% dplyr::select(-muttype_sub) %>% 
  dplyr::group_by(muttype) %>% dplyr::summarise_all(list(~sum(.))) %>% 
  tidyr::gather(key = "sample", value = "count", -.data$muttype) %>% 
  dplyr::mutate(sample = factor(sample, levels = unique(sample)))
# nr_indels <- indel_counts_main %>% dplyr::group_by(sample) %>% dplyr::summarise(nr_muts = sum(count))
# facet_labs_y_indels <- stringr::str_c(nr_muts$sample, " (n = ", 
#                                nr_muts$nr_muts, ")")
# names(facet_labs_y_indels) <- nr_muts$sample
indel_colors <- c("#FDBE6F", "#FF8001", "#B0DD8B", "#36A12E", "#FDCAB5", 
            "#FC8A6A", "#F14432", "#BC141A", "#D0E1F2", "#94C4DF", 
            "#4A98C9", "#1764AB", "#E2E2EF", "#B6B6D8", "#8683BD", 
            "#61409B")
indel_counts_main$label <- ifelse(indel_counts_main$count > 200, indel_counts_main$count, "")
indel_counts_main$genelabel <- gsub("-.*", "", indel_counts_main$sample)
indel_counts_main$genelabel <- factor(indel_counts_main$genelabel, levels = c("WT", "FANCCKO", "MSH2KO"))

indel_counts_main$samplelabel <- gsub(".*-", "", indel_counts_main$sample)
indel_counts_main$samplelabel <- gsub("A\n", "A-", indel_counts_main$samplelabel)
indel_counts_main$samplelabel <- gsub("\n", "", indel_counts_main$samplelabel)
indel_counts_main$samplelabel <- factor(indel_counts_main$samplelabel, levels = c("Subclone1","PTA-Raw","PTA-PTATO" ))

levels(indel_counts_main$muttype) <- gsub(levels(indel_counts_main$muttype), pattern = "deletion", replacement = "del")
levels(indel_counts_main$muttype) <- gsub(levels(indel_counts_main$muttype), pattern = "insertion", replacement = "ins")
levels(indel_counts_main$muttype) <- gsub(levels(indel_counts_main$muttype), pattern = "microhomology", replacement = "MH")

indel_counts_main_plot <- ggplot(indel_counts_main, aes(x = samplelabel, y = count, fill = muttype,  label = as.character(round(as.numeric(label),0)))) + 
  geom_bar(stat = "identity",  col = "black", width = 0.75, linewidth = 0.2) + 
  geom_text( size = 1.6, position = position_stack(vjust = 0.5)) +
  geom_text(aes(label = round(stat(y), 0), group = samplelabel), 
            stat = 'summary', fun = sum, vjust = -1, size = 1.6) +
  facet_grid(.~genelabel, space = "free", scales = "free", switch = "x") +
  labs(x = "", y = "Indels (#)", fill = "Indel type") + 
  scale_fill_manual(values = indel_colors) + 
  scale_y_continuous(expand = expansion(mult = c(0,0.1))) +   
  theme_classic(base_size = 6) + 
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), 
        strip.placement = "outside",
        legend.key.size = unit(0.2, 'cm'), 
        legend.margin = margin(t = 5, b = 5, l= 0, r =0, unit = "pt")) +
  guides(fill=guide_legend(ncol=2))

pdf(file = paste(Output_dir_FigureS6, "FigureS6D_IndelSpectra_abs.pdf", sep = ""), width = 3.8, height = 2, onefile = F)
plot_adjusted_strips(indel_counts_main_plot, strip_colors = AHH1_colors)
dev.off()


indel_counts_main_plot_relative <- ggplot(indel_counts_main, aes(x = samplelabel, y = count, fill = muttype,  label = as.character(round(as.numeric(label),0)))) + 
  geom_bar(stat = "identity",  col = "black", width = 0.75, linewidth = 0.2, position = position_fill()) + 
  facet_grid(.~genelabel, space = "free", scales = "free", switch = "x") +
  labs(x = "", y = "Indels (% of total)", fill = "Indel type") + 
  scale_fill_manual(values = indel_colors) + 
  scale_y_continuous(expand = expansion(mult = c(0,0.01)), labels = scales::percent) +   
  theme_classic(base_size = 6) + 
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), 
        strip.placement = "outside",
        legend.key.size = unit(0.2, 'cm'), 
        legend.margin = margin(t = 5, b = 5, l= 0, r =0, unit = "pt")) +
  guides(fill=guide_legend(ncol=2))

pdf(file = paste(Output_dir_FigureS6, "FigureS6E_IndelSpectra_rel.pdf", sep = ""), width = 3.8, height = 2, onefile = F)

plot_adjusted_strips(indel_counts_main_plot_relative, strip_colors = AHH1_colors)
dev.off()

# Determine which variants present in the preceding clone are also called by PTATO or SCAN2 in the PTA sample
Overview_CloneVariants <- data.frame()
ClonalVariants_SCAN2_overview <- data.frame()

for(Sample in Metadata$Sample[ Metadata$Type == "PTA"]){
 print(Sample) 
  Cell_line <- unique(Metadata$Cell_line[Metadata$Sample == Sample])
  
  Clone <- Metadata$Clone[Metadata$Sample == Sample]

  Samples <- Metadata$Sample[which(Metadata$Clone == Clone)]
  
  # Select the (sub)clone that preceded the PTA sample:
  Previous_Clone <- Metadata$Sample[Metadata$Days_after_clone == max(Metadata$Days_after_clone[Metadata$Sample %in% Samples & Metadata$Sample != Sample]) & Metadata$Sample %in% Samples]
  
  Variant_table <- variant_tables[[Cell_line]]
  SomaticVariants <- multisample_VCFs[[Cell_line]]
  
  sSNVs_sample <- SomaticVariants[isSNV(SomaticVariants)]
  sIndels_sample <- SomaticVariants[isIndel(SomaticVariants)]
  
  # Collect the SNVs that are CLONAL and PASS_QC in the (sub)clone preceding the PTA sample
  Clone_variants <- SomaticVariants[rownames(Variant_table)[Variant_table[,Previous_Clone] == 1 & !is.na(Variant_table[,Previous_Clone])]]
  Clone_variants_g <- rowRanges(Clone_variants)
  Clone_variants_g$VAF <- unlist(geno(Clone_variants)$VAF[,Previous_Clone])
  Clone_variants_g <- Clone_variants_g[Clone_variants_g$VAF > MinimalVAF] # The VAF cutoff for clonal variants in SMuRF may have been 0.15, which is too low 
  
  ## PASS_QC vs FAIL_QC
  ClonalVariants <- data.frame(Variant = names(Clone_variants_g))
  ClonalVariants$Type <-"Indel"
  ClonalVariants$Type[ClonalVariants$Variant %in% names(sSNVs_sample)] <- "SNV"
  
  # Variants with a NA value in the variant_table have a FAIL_QC label by SMuRF
  ClonalVariants$QC <- ifelse(is.na(Variant_table[ClonalVariants$Variant, Sample]) == T, "FAIL_QC", "PASS_QC")
  
  ### CALLABLELOCI: to determine if variants present in the clone have sufficient coverage in the PTA sample
  CallableLoci_g <- CallableLoci[[Sample]]
  
  olap_callable <- findOverlaps(Clone_variants_g, CallableLoci_g)
  Callable_variants <- names(Clone_variants_g)[queryHits(olap_callable)]
  
  ClonalVariants$CALLABLE <- ifelse(ClonalVariants$Variant %in% Callable_variants, "CALLABLE", "UNCALLABLE")
  
  grl_raw_sample <- sSNVs_raw[[Sample]]
  grl_PASS <- list()
  grl_FAIL <- list()
  grl_lowVAF <- list()
  # Remove all variants that are UNCALLABLE and/or with a PTAprobs < PTAprobsCutoff
  grl_PASS[[Sample]] <- grl_raw_sample[grl_raw_sample$FILTER == "PASS"]
  grl_FAIL[[Sample]] <- grl_raw_sample[grl_raw_sample$FILTER == "FAIL"]
  grl_lowVAF[[Sample]] <- grl_raw_sample[grl_raw_sample$FILTER == "FAIL_VAF"]
  
  # 
  sIndels_PASS <- sIndels_PTA_flag[[Sample]][sIndels_PTA_flag[[Sample]]$FILTER =="PASS"]
  sIndels_FAIL <- sIndels_PTA_flag[[Sample]][sIndels_PTA_flag[[Sample]]$FILTER !="PASS"]
  # sIndels_HOMOPOLYMER <- sIndels_PTA_flag[[Sample]][sIndels_PTA_flag[[Sample]]$FILTER =="HOMOPOLYMER"]
  # sIndels_REC_HOMOPOL <- sIndels_PTA_flag[[Sample]][sIndels_PTA_flag[[Sample]]$FILTER == "RECURRENT;HOMOPOLYMER"]
  
  
  # SNV
  ClonalVariants$PTATO <- NA
  ClonalVariants$PTATO[ClonalVariants$Variant %in% names(grl_PASS[[Sample]])] <- "PASS"
  ClonalVariants$PTATO[ClonalVariants$Variant %in% names(grl_FAIL[[Sample]])] <- "FAIL"
  ClonalVariants$PTATO[ClonalVariants$Variant %in% names(grl_lowVAF[[Sample]])] <- "FAIL_VAF"
  # Indels
  ClonalVariants$PTATO[ClonalVariants$Variant %in% names(sIndels_PASS)] <- "PASS"
  ClonalVariants$PTATO[ClonalVariants$Variant %in% names(sIndels_FAIL)] <- "FAIL"
  
  ClonalVariants$Conclusion <- ClonalVariants$PTATO
  
  # FAIL_QC and LOW_COV overwrite FAIL/PASS/RECURRENT etc
  ClonalVariants$Conclusion[ClonalVariants$QC != "PASS_QC"] <- "FAIL_QC"
  ClonalVariants$Conclusion[ClonalVariants$CALLABLE != "CALLABLE"] <- "LOW_COV"
  ClonalVariants$Sample <- Sample
  ClonalVariants$Caller <- "PTATO"
  
  ClonalVariants_SCAN2 <- ClonalVariants
  ClonalVariants_SCAN2$Conclusion <- NA
  ClonalVariants_SCAN2$Caller <- "SCAN2"
  ClonalVariants_SCAN2$Conclusion[ClonalVariants$QC != "PASS_QC"] <- "FAIL_QC"
  ClonalVariants_SCAN2$Conclusion[ClonalVariants$CALLABLE != "CALLABLE"] <- "LOW_COV"
  ClonalVariants_SCAN2$Conclusion[ClonalVariants_SCAN2$Variant %in% names(sSNVs_SCAN2[[Sample]])] <- "PASS"
  
  ClonalVariants_merged <- rbind(ClonalVariants, ClonalVariants_SCAN2)
  Overview_CloneVariants <- rbind(Overview_CloneVariants, ClonalVariants_merged)
  #ClonalVariants_SCAN2_overview <- rbind(ClonalVariants_SCAN2_overview, ClonalVariants_SCAN2)
}

# SCAN2 didnt run for chr17, exclude
Overview_CloneVariants <- Overview_CloneVariants[startsWith(prefix = "17:", Overview_CloneVariants$Variant) == FALSE,]

Overview_CloneVariants$ID <- Overview_CloneVariants$Sample
Overview_CloneVariants$ID[Overview_CloneVariants$Sample == "PMCAHH1-WT-C6SC1"] <- "WT-\nPTA1"
Overview_CloneVariants$ID[Overview_CloneVariants$Sample == "PMCAHH1-WT-C19SC1"] <- "WT-\nPTA2"
Overview_CloneVariants$ID[Overview_CloneVariants$Sample == "PMCAHH1-FANCCKO-C02B03SC03E05-PTAP1D7"] <- "FANCCKO-\nPTA1"
Overview_CloneVariants$ID[Overview_CloneVariants$Sample == "PMCAHH1-MSH2KO-C27E06SC51B06-PTAP1E7"] <- "MSH2KO-\nPTA1"
Overview_CloneVariants$ID <- factor(Overview_CloneVariants$ID, levels = c("WT-\nPTA1", "WT-\nPTA2", "FANCCKO-\nPTA1", "MSH2KO-\nPTA1"))
# SNV_Counts <- Overview_CloneVariants[Overview_CloneVariants$Type == "SNV",] %>% group_by(ID, Conclusion, Caller) %>%
#   summarise(n = n()) %>%
#   mutate(freq = n / sum(n))
# SNV_Counts_Total <-  Overview_CloneVariants[Overview_CloneVariants$Type == "SNV",] %>% group_by(ID, Type, Caller) %>%
#   summarise(n = n()) %>%
#   mutate(freq = n / sum(n))

SNV_Counts_PTATO <- Overview_CloneVariants[Overview_CloneVariants$Type == "SNV" & Overview_CloneVariants$Caller == "PTATO",] %>% group_by(ID, Conclusion) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
SNV_Counts_PTATO$Caller <- "PTATO"
SNV_Counts_Total_PTATO <-  Overview_CloneVariants[Overview_CloneVariants$Type == "SNV" & Overview_CloneVariants$Caller == "PTATO",] %>% group_by(ID, Type) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
SNV_Counts_Total_PTATO$Caller <- "PTATO"

SNV_Counts_SCAN2 <- Overview_CloneVariants[Overview_CloneVariants$Type == "SNV" & Overview_CloneVariants$Caller == "SCAN2",] %>% group_by(ID, Conclusion) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
SNV_Counts_SCAN2$Caller <- "SCAN2"

SNV_Counts_Total_SCAN2  <- Overview_CloneVariants[Overview_CloneVariants$Type == "SNV" & Overview_CloneVariants$Caller == "SCAN2",] %>% group_by(ID, Type) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
SNV_Counts_Total_SCAN2$Caller <- "SCAN2"

SNV_Counts <- rbind(SNV_Counts_PTATO, SNV_Counts_SCAN2)
SNV_Counts_Total <- rbind(SNV_Counts_Total_PTATO, SNV_Counts_Total_SCAN2)


SNV_Counts$Label <- round(SNV_Counts$freq, 2)
SNV_Counts$Label <- paste(SNV_Counts$Label*100 , "%", sep = "")
SNV_Counts$Label[SNV_Counts$freq < 0.08] <- ""
SNV_Counts$Conclusion[is.na(SNV_Counts$Conclusion)] <- "ABSENT"
SNV_Counts$Conclusion <- factor(SNV_Counts$Conclusion, levels = c("ABSENT", "LOW_COV", "FAIL_QC","FAIL","PASS"))

ClassificationColors <- c("darkgray", "#7570B3","#E6AB02", "#D55E00", "#009E73" )
names(ClassificationColors) <- c("ABSENT", "LOW_COV", "FAIL_QC","FAIL","PASS")

SNV_Counts_plot <- ggplot(SNV_Counts, aes(x = Caller, y = freq, fill = Conclusion)) +
  geom_bar(stat = "identity", width = 0.85, col = "black", linewidth = 0.2) +
  facet_grid(.~ID, space = "free", scales = "free", switch = "x") +
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 1.8, col = "white") +
  geom_text(data = SNV_Counts_Total, aes(label = n, y = 1.05, x = Caller), inherit.aes = F, size = 1.8) +
  theme_classic(base_size = 6) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.1), labels = scales::percent) + 
  scale_fill_manual(values = ClassificationColors) +
  #scale_fill_manual(values = c("#D55E00","#009E73","#999999")) +
  labs(fill = "Filter", x = "Sample", y = "Shared base substitutions\n(%)") +
  theme(axis.text.x = element_text(angle = 45,  hjust=1),
        strip.placement = "outside",
         strip.background = element_rect(color = "white"),
        # strip.text.x = element_blank(),
        legend.key.size = unit(0.25, 'cm')) 

#pdf(file =paste(Output_dir_Figure2, "Fig2f_Shared_BaseSubs.pdf",sep = ""), width = 2.75, height = 1.75, onefile = F)

plot_adjusted_strips(SNV_Counts_plot, strip_colors = c(AHH1_colors[c(1,1,2,3)]))
#dev.off()

# Average the 4 cell lines
SNV_Counts_mean <- SNV_Counts  %>% group_by(Caller, Conclusion) %>%
  summarise(across(freq, list(mean = mean, sd = sd)))
SNV_Counts_mean$Label <- paste(round(SNV_Counts_mean$freq_mean*100, 0), "%", sep = "")

SNV_Counts_Total_mean <- SNV_Counts_Total  %>% group_by(Caller) %>%
  summarise(across(n, list(sum = sum)))

ggplot(SNV_Counts_mean, aes(x = Caller, y = freq_mean, fill = Conclusion)) +
  geom_bar(stat = "identity", width = 0.85, col = "black", linewidth = 0.2) +
  #geom_errorbar( aes(ymin=freq_mean-freq_sd, ymax=freq_mean+freq_sd), width=.2,position = "identity")  +
  
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 1.8, col = "white") +
  geom_text(data = SNV_Counts_Total_mean, aes(label = n_sum, y = 1.05, x = Caller), inherit.aes = F, size = 1.8) +
  theme_classic(base_size = 6) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.1), labels = scales::percent) + 
  scale_fill_manual(values = ClassificationColors) +
  #scale_fill_manual(values = c("#D55E00","#009E73","#999999")) +
  labs(fill = "Filter", x = "Sample", y = "Shared base substitutions\n(%)") +
  theme(axis.text.x = element_text(angle = 45,  hjust=1),
        strip.placement = "outside",
        strip.background = element_rect(color = "white"),
        # strip.text.x = element_blank(),
        legend.key.size = unit(0.25, 'cm')) 

ggsave(paste(Output_dir_Figure2, "Fig_2I_PMCAHH1_SharedBaseSubs.pdf", sep = ""), width = 4.5, height = 3.75, units = "cm")


CallableVariants <- Overview_CloneVariants$Variant[Overview_CloneVariants$Conclusion %in% c("PASS", "FAIL") & Overview_CloneVariants$Type == "SNV"  & Overview_CloneVariants$Caller == "PTATO"] 

SNV_Counts_Callable_PTATO <- Overview_CloneVariants[Overview_CloneVariants$Variant %in% CallableVariants & Overview_CloneVariants$Caller == "PTATO",] %>% group_by(ID, Conclusion) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
SNV_Counts_Callable_PTATO$Caller <- "PTATO"
SNV_Counts_Callable_Total_PTATO <- Overview_CloneVariants[Overview_CloneVariants$Variant %in% CallableVariants & Overview_CloneVariants$Caller == "PTATO",] %>% group_by(ID) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
SNV_Counts_Callable_Total_PTATO$Caller <- "PTATO"

SNV_Counts_Callable_SCAN2 <- Overview_CloneVariants[Overview_CloneVariants$Variant %in% CallableVariants & Overview_CloneVariants$Caller == "SCAN2",] %>% group_by(ID, Conclusion) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
SNV_Counts_Callable_SCAN2$Caller <- "SCAN2"

SNV_Counts_Callable_Total_SCAN2 <- Overview_CloneVariants[Overview_CloneVariants$Variant %in% CallableVariants & Overview_CloneVariants$Caller == "SCAN2",] %>% group_by(ID) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
SNV_Counts_Callable_Total_SCAN2$Caller <- "SCAN2"

SNV_Counts_Callable <- rbind(SNV_Counts_Callable_PTATO, SNV_Counts_Callable_SCAN2)
SNV_Counts_Callable_Total <- rbind(SNV_Counts_Callable_Total_PTATO, SNV_Counts_Callable_Total_SCAN2)

SNV_Counts_Callable$Label <- round(SNV_Counts_Callable$freq, 2)
SNV_Counts_Callable$Label  <- paste(SNV_Counts_Callable$Label * 100, "%", sep ="")
SNV_Counts_Callable$Label[SNV_Counts_Callable$freq < 0.05] <- ""

SNV_Counts_Callable$Conclusion[is.na(SNV_Counts_Callable$Conclusion)] <- "ABSENT"
SNV_Counts_Callable$Conclusion <- factor(SNV_Counts_Callable$Conclusion, levels = c("ABSENT", "LOW_COV", "FAIL_QC","FAIL","PASS"))


SNV_Counts_Callable_plot <- ggplot(SNV_Counts_Callable, aes(x = Caller, y = freq, fill = Conclusion)) + 
  geom_bar(stat = "identity", width = 0.85, col = "black", size = 0.2) +
  facet_grid(.~ID, space = "free", scales = "free", switch= "x") +
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 1.8, col = "white") +
  geom_text(data = SNV_Counts_Callable_Total, aes(label = n, y = 1.05, x = Caller), inherit.aes = F, size = 1.8) +
  theme_classic(base_size = 6) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.1) , labels = scales::percent) + 
  scale_fill_manual(values = ClassificationColors, drop = F) +
  #scale_fill_manual(values = c("#D55E00","#009E73","#999999")) +
  labs(fill = "Filter", x = "Sample", y = "Shared base substitutions\nin callable regions (%)") +
  theme(axis.text.x = element_text(angle = 45,  hjust=1),
        strip.placement = "outside",
        strip.background = element_rect(color = "white"),
        # strip.text.x = element_blank(),
        legend.key.size = unit(0.25, 'cm')) 

#pdf(file =paste(Output_dir_Figure2, "Fig2f_bottom_Shared_BaseSubs_Callable.pdf",sep = ""), width = 2.75, height = 1.75, onefile = F)

plot_adjusted_strips(SNV_Counts_Callable_plot, strip_colors = c(AHH1_colors[c(1,1,2,3)]))
#dev.off()

SNV_Counts_Callable_mean <- SNV_Counts_Callable  %>% group_by(Caller, Conclusion) %>%
  summarise(across(freq, list(mean = mean, sd = sd)))
SNV_Counts_Callable_mean$Label <- paste(round(SNV_Counts_Callable_mean$freq_mean*100, 0), "%", sep = "")
SNV_Counts_Callable_mean$Conclusion <- factor(SNV_Counts_Callable_mean$Conclusion, levels = c("ABSENT", "LOW_COV", "FAIL_QC","FAIL","PASS"))
SNV_Counts_Callable_Total_mean <- SNV_Counts_Callable_Total  %>% group_by(Caller) %>%
  summarise(across(n, list(sum = sum)))

ggplot(SNV_Counts_Callable_mean, aes(x = Caller, y = freq_mean, fill = Conclusion)) +
  geom_bar(stat = "identity", width = 0.85, col = "black", linewidth = 0.2) +
  #geom_errorbar( aes(ymin=freq_mean-freq_sd, ymax=freq_mean+freq_sd), width=.2,position = "identity")  +
  
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 1.8, col = "white") +
  geom_text(data = SNV_Counts_Callable_Total_mean, aes(label = n_sum, y = 1.05, x = Caller), inherit.aes = F, size = 1.8) +
  theme_classic(base_size = 6) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.1), labels = scales::percent) + 
  scale_fill_manual(values = ClassificationColors, drop = F) +
  #scale_fill_manual(values = c("#D55E00","#009E73","#999999")) +
  labs(fill = "Filter", x = "Sample", y = "Shared base substitutions\nin callable regions (%)") +
  theme(axis.text.x = element_text(angle = 45,  hjust=1),
        strip.placement = "outside",
        strip.background = element_rect(color = "white"),
        # strip.text.x = element_blank(),
        legend.key.size = unit(0.25, 'cm')) 

ggsave(paste(Output_dir_Figure2, "Fig_2J_PMCAHH1_SharedBaseSubs_Callable.pdf", sep = ""), width = 4.5, height = 3.75, units = "cm")



### Indels
Indel_Counts_PTATO <- Overview_CloneVariants[Overview_CloneVariants$Type == "Indel" & Overview_CloneVariants$Caller == "PTATO",] %>% group_by(ID, Conclusion) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
Indel_Counts_PTATO$Caller <- "PTATO"
Indel_Counts_Total_PTATO <-  Overview_CloneVariants[Overview_CloneVariants$Type == "Indel" & Overview_CloneVariants$Caller == "PTATO",] %>% group_by(ID, Type) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
Indel_Counts_Total_PTATO$Caller <- "PTATO"


Indel_Counts_PTATO$Label <- round(Indel_Counts_PTATO$freq, 2)
Indel_Counts_PTATO$Label <- paste(Indel_Counts_PTATO$Label*100 , "%", sep = "")
Indel_Counts_PTATO$Label[Indel_Counts_PTATO$freq < 0.08] <- ""
Indel_Counts_PTATO$Conclusion[is.na(Indel_Counts_PTATO$Conclusion)] <- "ABSENT"
#Indel_Counts_PTATO$Conclusion <- factor(Indel_Counts_PTATO$Conclusion, levels = c("ABSENT", "LOW_COV", "FAIL_QC","FAIL","PASS"))

ClassificationColorsIndels <- c("darkgray", "#7570B3","#E6AB02", "#D55E00", "#009E73" )
names(ClassificationColorsIndels) <- c("ABSENT", "LOW_COV", "FAIL_QC","FAIL","PASS")


Indel_Counts_plot <- ggplot(Indel_Counts_PTATO, aes(x = Caller, y = freq, fill = Conclusion)) +
  geom_bar(stat = "identity", col = "black", size = 0.2) +
  facet_grid(.~ID, space = "free", scales = "free", switch = "x") +
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 1.8, col = "white") +
  geom_text(data = Indel_Counts_Total_PTATO, aes(label = n, y = 1.05, x = Caller), inherit.aes = F, size = 1.8) +
  theme_classic(base_size = 6) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.1), labels = scales::percent) + 
  scale_fill_manual(values = ClassificationColorsIndels) +
  labs(fill = "Filter", x = "Caller", y = "Shared indels\n(%)") +
  theme(axis.text.x = element_blank(),
        strip.background = element_rect(color = "white"),
        strip.text.x = element_text(angle = 90),
        legend.key.size = unit(0.25, 'cm'),
        strip.placement = "outside") 

#pdf(file =paste(Output_dir_Figure2, "Fig2g_Shared_Indels.pdf",sep = ""), width = 2, height = 1.75, onefile = F)
plot_adjusted_strips(Indel_Counts_plot, strip_colors = c(AHH1_colors[c(1,1,2,3)]))
#dev.off()

# Average the 4 cell lines
Indel_Counts_mean <- Indel_Counts_PTATO  %>% group_by(Caller, Conclusion) %>%
  summarise(across(freq, list(mean = mean, sd = sd)))
Indel_Counts_mean$Label <- paste(round(Indel_Counts_mean$freq_mean*100, 0), "%", sep = "")
Indel_Counts_Total_PTATO_mean <- Indel_Counts_Total_PTATO  %>% group_by(Caller) %>%
  summarise(across(n, list(sum = sum)))

Indel_Counts_mean_plot <- ggplot(Indel_Counts_mean, aes(x = Caller, y = freq_mean, fill = Conclusion)) +
  geom_bar(stat = "identity", col = "black", size = 0.2) +
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 1.8, col = "white") +
  geom_text(data = Indel_Counts_Total_PTATO_mean, aes(label = n_sum, y = 1.05, x = Caller), inherit.aes = F, size = 1.8) +
  theme_classic(base_size = 6) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.1), labels = scales::percent) + 
  scale_fill_manual(values = ClassificationColorsIndels) +
  labs(fill = "Filter", x = "Caller", y = "Shared indels\n(%)") +
  theme(legend.key.size = unit(0.25, 'cm'), axis.text.x = element_text(angle = 45,  hjust=1))
Indel_Counts_mean_plot

ggsave(paste(Output_dir_Figure2, "Fig_2K_PMCAHH1_SharedIndels.pdf", sep = ""), width = 4, height = 3.5, units = "cm")

  # Select only the indels that are in callable regions, with sufficient quality
CallableIndels<- Overview_CloneVariants$Variant[Overview_CloneVariants$Conclusion %in% c("PASS", "FAIL") & Overview_CloneVariants$Type == "Indel"  & Overview_CloneVariants$Caller == "PTATO"] 

Indel_Counts_Callable_PTATO <- Overview_CloneVariants[Overview_CloneVariants$Variant %in% CallableIndels & Overview_CloneVariants$Caller == "PTATO",] %>% group_by(ID, Conclusion) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
Indel_Counts_Callable_PTATO$Caller <- "PTATO"
Indel_Counts_Callable_Total_PTATO <- Overview_CloneVariants[Overview_CloneVariants$Variant %in% CallableIndels & Overview_CloneVariants$Caller == "PTATO",] %>% group_by(ID) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
Indel_Counts_Callable_Total_PTATO$Caller <- "PTATO"

Indel_Counts_Callable_PTATO$Label <- round(Indel_Counts_Callable_PTATO$freq, 2)
Indel_Counts_Callable_PTATO$Label  <- paste(Indel_Counts_Callable_PTATO$Label * 100, "%", sep ="")
Indel_Counts_Callable_PTATO$Label[Indel_Counts_Callable_PTATO$freq < 0.05] <- ""

Indel_Counts_Callable_PTATO$Conclusion[is.na(Indel_Counts_Callable_PTATO$Conclusion)] <- "ABSENT"
Indel_Counts_Callable_PTATO$Conclusion <- factor(Indel_Counts_Callable_PTATO$Conclusion, levels = c("ABSENT", "LOW_COV", "FAIL_QC","FAIL","PASS"))

Indels_Counts_Callable_plot <- ggplot(Indel_Counts_Callable_PTATO, aes(x = Caller, y = freq, fill = Conclusion)) + 
  geom_bar(stat = "identity", col = "black", size = 0.2) +
  facet_grid(.~ID, space = "free", scales = "free", switch= "x") +
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 1.8, col = "white") +
  geom_text(data = Indel_Counts_Callable_Total_PTATO, aes(label = n, y = 1.05, x = Caller), inherit.aes = F, size = 1.8) +
  theme_classic(base_size = 6) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.1) , labels = scales::percent) + 
  scale_fill_manual(values = ClassificationColorsIndels, drop = F) +
  #scale_fill_manual(values = c("#D55E00","#009E73","#999999")) +
  labs(fill = "Filter", x = "Sample", y = "Shared callable indels\nin callable regions (%)") +
  theme(axis.text.x = element_blank(),
        strip.background = element_rect(color = "white"),
        strip.text.x = element_text(angle = 90),
        legend.key.size = unit(0.25, 'cm'),
        strip.placement = "outside") 

#pdf(file =paste(Output_dir_Figure2, "Fig2g_Shared_Indels_Callable.pdf",sep = ""), width = 2, height = 1.75, onefile = F)
plot_adjusted_strips(Indels_Counts_Callable_plot, strip_colors = c(AHH1_colors[c(1,1,2,3)]))
#dev.off()

Indel_Counts_Callable_mean <- Indel_Counts_Callable_PTATO  %>% group_by(Caller, Conclusion) %>%
  summarise(across(freq, list(mean = mean, sd = sd)))
Indel_Counts_Callable_mean$Label <- paste(round(Indel_Counts_Callable_mean$freq_mean*100, 0), "%", sep = "")
Indel_Counts_Callable_mean$Conclusion <- factor(Indel_Counts_Callable_mean$Conclusion, levels = c("ABSENT", "LOW_COV", "FAIL_QC","FAIL","PASS"))
Indel_Counts_Callable_Total_mean <- Indel_Counts_Callable_Total_PTATO  %>% group_by(Caller) %>%
  summarise(across(n, list(sum = sum)))


ggplot(Indel_Counts_Callable_mean, aes(x = Caller, y = freq_mean, fill = Conclusion)) +
  geom_bar(stat = "identity", width = 0.85, col = "black", linewidth = 0.2) +
  #geom_errorbar( aes(ymin=freq_mean-freq_sd, ymax=freq_mean+freq_sd), width=.2,position = "identity")  +
  
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 1.8, col = "white") +
  geom_text(data = Indel_Counts_Callable_Total_mean, aes(label = n_sum, y = 1.05, x = Caller), inherit.aes = F, size = 1.8) +
  theme_classic(base_size = 6) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.1), labels = scales::percent) + 
  scale_fill_manual(values = ClassificationColors, drop = F) +
  #scale_fill_manual(values = c("#D55E00","#009E73","#999999")) +
  labs(fill = "Filter", x = "Sample", y = "Shared indels\nin callable regions (%)") +
  theme(axis.text.x = element_text(angle = 45,  hjust=1),
        legend.key.size = unit(0.25, 'cm')) 

ggsave(paste(Output_dir_Figure2, "Fig_2L_PMCAHH1_SharedIndels_Callable.pdf", sep = ""), width = 4, height = 3.5, units = "cm")

