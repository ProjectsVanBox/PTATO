### Figure 1, parts of S1 and S2
### First part contains code to plot base substitutions, second part contains code to plot indels

library(randomForest)
library(BSgenome.Hsapiens.UCSC.hg38)
library(gridExtra)
# library(cowplot)
library(ggpubr)
library(tibble)
library(MutationalPatterns)
library(reshape2)
library(VennDiagram)
library(RColorBrewer)
library(ggplot2)
library(ggh4x)
library(BSgenome)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

options(scipen = 999)

Input_dir <-  "/path/to/mendeleydata/"

Resources_dir <- paste(Input_dir, "Resources/", sep = "")
Scripts_dir <- "/path/to/scriptsgithub/"

Output_dir_Fig1 <-  "/path/to/dir/Figure1/"
Output_dir_FigS1 <- "/path/to/dir/FigureS1/"
Output_dir_FigS2 <-  "/path/to/dir/FigureS2/"

if(dir.exists(Output_dir_Fig1) == F){dir.create(Output_dir_Fig1)}
if(dir.exists(Output_dir_FigS1) == F){dir.create(Output_dir_FigS1)}
if(dir.exists(Output_dir_FigS2) == F){dir.create(Output_dir_FigS2)}

source(paste(Scripts_dir, "manual_training_performance_functions.R", sep = ""))
source(paste(Scripts_dir, "GeneralFunctions.R", sep = ""))
source(paste(Scripts_dir, "MutationalPatterns_Addon.R", sep = ""))

## Figure 1b top panel: Venn diagram training set
# This is just a visual representation of the training set, so no realistic representation (eg the numbers do not necessarily correspond to the real numbers in the training set)
grid.newpage()
pdf(file = paste(Output_dir_Fig1,"Fig1b_Venn_TPs.pdf", sep = ""), width = 0.6, height = 0.6, pointsize = 5)
draw.pairwise.venn(area1=900, area2=900,cross.area=756,
                   category=c("Bulk","PTA"),fill=c("#E69F00","#0072B2"),
                   fontfamily = "sans",
                   cat.fontfamily = "sans", lwd = 1, cex = 0.8)
dev.off()

grid.newpage()
pdf(file = paste(Output_dir_Fig1, "Fig1b_Venn_FPs.pdf", sep = ""), width = 0.6, height = 0.6, pointsize = 5)
draw.triple.venn(area1=2000, area2=2000, area3=2000, 
                 n12=250, n23=250, n13=250, n123=200, 
                 category=c("Cell 1","Cell 2","Cell 3"),
                 col="black",fill=c("#D55E00","#E69F00","#F0E442"),
                 fontfamily = "sans",
                 cat.fontfamily = "sans", lwd = 1, cex = 0.8)
dev.off()


### Figure 1b: Number of variants used in the training set per sample


false_pos_fnames = list.files(paste(Input_dir, "Training/intermediate/snvs/rf/", sep = ""), pattern = "lowwalker|other", recursive = TRUE, full.names = TRUE)
true_pos_fnames = list.files(paste(Input_dir, "Training/intermediate/snvs/rf/", sep = ""), pattern = "shared", recursive = TRUE, full.names = TRUE)
nopta_label = "noPTA"
pta_label = "PTA"

# Read and combine training data
#train.data.1.l <-lapply(strsplit(shared_fnames,",")[[1]], read_train_data, nopta_label)
#train.data.2.l <-lapply(strsplit(lowwalker_fnames,",")[[1]], read_train_data, pta_label)
train.data.1.l <-lapply(true_pos_fnames, read_train_data, nopta_label)
train.data.2.l <-lapply(false_pos_fnames, read_train_data, pta_label)
train.data <- do.call(rbind, c(train.data.1.l, train.data.2.l))

train.data = dplyr::filter(train.data, DONOR_ID != "PMCAHH1-MSH2KO")

# Clean up data
train.data <- train.data[!(duplicated(train.data[,c("CHROM","START","END","val")])),]
train.data <- train.data[!(duplicated(train.data[,c("CHROM","START","END")]) | duplicated(train.data[,c("CHROM","START","END")], fromLast = T)), ]

train.data.1 <- na.omit(train.data)
train.data.1 <- train.data.1[is.finite(train.data.1$p_binom),]
train.data.1 = subset_classes_equal(train.data.1)

nr_muts_donor = train.data.1[,c("DONOR_ID", "val", "POS0")] %>% 
  table() %>% 
  as.data.frame()
levels(nr_muts_donor$DONOR_ID)[levels(nr_muts_donor$DONOR_ID) == "PMCAHH1-FANCCKO"] <- "PMCAHH1-\nFANCCKO"


nr_muts_donor_total = train.data.1[,c("DONOR_ID", "val", "POS0")] %>% 
  table() %>% 
  as.data.frame()

ggplot(nr_muts_donor, aes(x = DONOR_ID, y = Freq, fill = POS0)) + 
  geom_bar(stat = "identity", width = 0.8) + facet_grid(.~val) + scale_fill_manual(values = COLORS6) +
  theme_classic(base_size = 6) + scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.25, 'cm'),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  labs(x = "Donor ID", y = "Base substitutions\nin training set (#)", fill = "Mutation type")
ggsave(paste(Output_dir_Fig1, "Fig1b_Training_nr_muts_door.pdf", sep = ""), width = 7.5, height = 3.5, units = "cm")


## Figure 1c
# Random forest feature importance
RF <- readRDS(paste(Input_dir, "Training/randomforest_v1.0.0.rds", sep = ""))
Feature_Importance = importance(RF$all) %>% 
  as.data.frame() %>% 
  rownames_to_column()

Feature_Importance$rowname[Feature_Importance$rowname == "p_binom"] <- "ALLELIC IMBALANCE"
Feature_Importance$rowname[Feature_Importance$rowname == "POS0"] <- "MUTATION TYPE"
Feature_Importance$rowname[Feature_Importance$rowname == "REPLISEQ"] <- "REPLICATION TIMING"

Feature_Importance$rowname <- gsub(x = Feature_Importance$rowname, pattern = "POSm", replacement = "POSITION -")
Feature_Importance$rowname <- gsub(x = Feature_Importance$rowname, pattern = "POSp", replacement = "POSITION +")
Feature_Importance_top10 <- Feature_Importance[order(Feature_Importance$MeanDecreaseGini, decreasing = T),][1:10,]

ggplot(Feature_Importance_top10, aes(x = reorder(rowname, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "#1F78B4") +
  coord_flip() +
  labs(x = "Features (top 10)", y = "Feature importance\n(Mean decrease Gini)") +
  scale_y_continuous(expand = c(0,0), limits = c(0,135)) +
  theme_bw(base_size = 6) +
  theme(axis.ticks.x = element_line(colour = "black", size = 0.25),
        axis.ticks.y = element_line(colour = "black", size = 0.25),
        axis.ticks.length=unit(0.05, "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
ggsave(paste(Output_dir_Fig1, "Fig1c_FeatureImportance.pdf", sep = ""), width = 5, height = 3.5, units = "cm")

## Figure 1d
# Precision-Recall
cutoff_accuracies <- readRDS(paste(Input_dir, "Training/cutoff_accuracies.rds", sep = ""))
auc <- calculate_auc(cutoff_accuracies)
roc_fig <- ggplot(cutoff_accuracies, aes(y = precision, x = recall)) +
  geom_line() +
  geom_point(col = "#1F78B4", size = 0.4, alpha = 0.8) +
  geom_area(alpha = 0.1, position = 'identity') +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) +
  labs(x = "Recall", y = "Precision") +
  annotate("text", label = paste0("AUC: ", round(auc, 3)), x = 0.25, y = 0.1, size = 1.6) +
  theme_bw(base_size = 6) +
  theme(axis.ticks.x = element_line(colour = "black", size = 0.25),
        axis.ticks.y = element_line(colour = "black", size = 0.25),
        axis.ticks.length=unit(0.05, "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
roc_fig
ggsave(paste(Output_dir_Fig1, "Fig1d_Training_AUC.pdf", sep = ""), width = 4, height = 3, units = "cm")

## Figure 1e
PTA_signature <- read.table(paste(Resources_dir, "/PTA_Artefact_Signature.txt", sep = ""), header=T)

value_mutmat = readRDS(paste(Input_dir, "Training/value_mutpatterns.rds",sep  = ""))[[1]]
colnames(value_mutmat) <- c("FP", "TP")
predicted_mutmat = readRDS(paste(Input_dir, "Training/predicted_mutpatterns.rds",sep  = ""))[[1]]
colnames(predicted_mutmat) <- c("FP", "TP")

dev.off()
pdf(paste(Output_dir_Fig1, "Fig1e_Training_96profile_input.pdf", sep = ""), width = 3, height = 1.5, onefile = F)
plot_96_profile_small(value_mutmat[,c("TP", "FP")], condensed = TRUE, ymax = 0.08)
dev.off()

pdf(paste(Output_dir_Fig1, "Fig1e_Training_96profile_output.pdf", sep = ""), width = 3, height = 1.5, onefile = F)
plot_96_profile_small(predicted_mutmat[,c("TP", "FP")], condensed = TRUE, ymax = 0.08)
dev.off()


### Figure 1f: Cosine similarities RF input and output
predicted_mutmat2 <- cbind(predicted_mutmat, PTA_signature)
colnames(predicted_mutmat2) <- c("FP", "TP", "PTA sig")

colnames(predicted_mutmat) <- c("FP", "TP")
colnames(value_mutmat) <- c("FP", "TP")
cos_sim_before_after <- melt(cos_sim_matrix(value_mutmat, predicted_mutmat2))

colnames(cos_sim_before_after) <- c("Input", "Output", "CosineSimilarity")
cos_sim_before_after$Input <- factor(cos_sim_before_after$Input, levels = c("FP", "TP"))
cos_sim_before_after$Output <- factor(cos_sim_before_after$Output, levels = c( "TP", "FP", "PTA sig"))
cos_sim_before_after$Type <- ifelse(cos_sim_before_after$Output == "PTA sig", "PTA", "RF")
cos_sim_before_after$Type <- factor(cos_sim_before_after$Type, levels = c("RF","PTA sig"))

ggplot(cos_sim_before_after, aes(x = Output, y = Input, fill = CosineSimilarity)) + geom_tile() +
  geom_text(aes(label = round(CosineSimilarity, 2)), size = 2) + 
  facet_grid(.~Type, scales= "free", space = "free") +
  scale_x_discrete(position = "top", expand = c(0,0))  +
  scale_y_discrete(expand = c(0,0))  +
  scale_fill_gradient(
    low = "cornsilk",
    high = "#1F78B4",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill") +
  labs(fill = "Cosine\nSimilarity", x = "Random forest output", y = "Random forest input") +
  theme_classic(base_size = 6) +
  theme(        strip.background = element_blank(),
                strip.text.x = element_blank(),
                legend.key.height= unit(0.3, 'cm'),
                legend.key.width= unit(0.3, 'cm'),
                legend.position = "bottom")
ggsave(paste(Output_dir_Fig1, "Fig1f_Training_InputOutput_CosineSim.pdf", sep =""), width = 3.5, height = 4, units = "cm")

### Figure 1D
snv_conf <- read.table(paste(Input_dir, "Training/randomforest_v1.0.0_confusion.txt", sep = ""), check.names = F, header = T)
snv_conf <- snv_conf[,1:2]

df1 <- data.frame(real = c("noPTA\n(Pos)",
                    "PTA\n(Neg)", "noPTA\n(Pos)", "PTA\n(Neg)"), 
           predicted = c("noPTA\n(Pos)",
                        "PTA\n(Neg)", "PTA\n(Neg)", "noPTA\n(Pos)"),
           values = c(snv_conf["noPTA", "noPTA"], 
                      snv_conf["PTA", "PTA"], 
                      snv_conf["noPTA", "PTA"], 
                      snv_conf["PTA", "noPTA"]),
           label = c("TP", "TN", "FN", "FP"))
df1$rel <- df1$values / sum(df1$values)

df2 <- data.frame(real = c("Accuracy", "Accuracy"), 
                  predicted = c("noPTA\n(Pos)", "PTA\n(Neg)"),
                  values = c( df1$values[df1$label == "FP"] / sum(df1$values[df1$predicted == "noPTA\n(Pos)"] ),
                              df1$values[df1$label == "FN"] / sum(df1$values[df1$predicted == "PTA\n(Neg)"])),
                  label = c("FDR", "FOR"))
df2$rel <- df2$values

df3 <- data.frame(real = c("noPTA\n(Pos)", "PTA\n(Neg)", "Accuracy"), 
                  predicted = c("Negatives", "Negatives", "Negatives"),
                  values = c( df1$values[df1$label == "TP"] / sum(df1$values[df1$real == "noPTA\n(Pos)"] ),
                              df1$values[df1$label == "TN"] / sum(df1$values[df1$real == "PTA\n(Neg)"]),
                              (df1$values[df1$label == "TP"]+df1$values[df1$label == "TN"]) / (sum(df1$values[df1$real == "noPTA\n(Pos)"]) + sum(df1$values[df1$real == "PTA\n(Neg)"]))),
                  label = c("TPR", "TNR", "ACCU"))
df3$rel <- df3$values



df_merged <- rbind(df1, df2, df3)

df_merged$label2 <- paste(df_merged$label, "\n", round(df_merged$values,2), sep = "")
df_merged$real <- factor(df_merged$real, levels = c( "Accuracy","PTA\n(Neg)", "noPTA\n(Pos)"))
df_merged$predicted <- factor(df_merged$predicted, levels = c("noPTA\n(Pos)", "PTA\n(Neg)", "Negatives"))

ggplot(df_merged, aes(x = predicted, y = real, fill = rel, label = label2)) + 
  geom_tile(col = "white", linewidth = 0.5) + 
  geom_text(size = 1.8) +
  theme_minimal(base_size = 6 ) +
  labs(fill = "Ratio", y = "Real value", x = "Predicted value") +
  scale_fill_gradient(
    low = "cornsilk",
    high = "#1F78B4",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill") +
  scale_x_discrete(position = "top") + theme(legend.key.size= unit(0.3, 'cm'))

ggsave(paste(Output_dir_Fig1, "Fig1d_ConfusionMatrix.pdf", sep =""), width = 4.5, height = 3, units = "cm")


#### Indels
Metadata_file <- paste(Input_dir, "Table_S1.txt", sep = "")
Metadata <- read.delim(Metadata_file, header = T)
Metadata <- Metadata[,-which(colnames(Metadata) %in% c("Label2","Mutation","CallableLoci_autosomal"))]
Metadata <- Metadata[Metadata$Training == "Yes",]

indels_raw <- list()
Chromosomes <- c(1:22)

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

# Read the SMuRF VCFs containing the somatic indel
for(Individual in unique(Metadata$Individual)){
  print(Individual)
  VCFs_individual <- list.files(paste(Input_dir, Individual, "/PTATO/intermediate/short_variants/somatic_vcfs/", Individual, "/", sep = ""), full.names = T, pattern = ".filtered.vcf.gz$")
  for(Sample in Metadata$Sample[Metadata$Individual == Individual]){
    print(Sample)
    VCF_sample <- VCFs_individual[grep(Sample, x = VCFs_individual)]
    print(VCF_sample)
    
    SomaticVariants_Sample <- readVcf(VCF_sample, "hg38")
    SomaticVariants_Sample <- SomaticVariants_Sample[as.vector(seqnames(SomaticVariants_Sample)) %in% Chromosomes,] 
    
    Indels_Sample <- SomaticVariants_Sample[isIndel(SomaticVariants_Sample)]
    Indel_rr <- VCF_to_GR(Indels_Sample, chromosomes = Chromosomes)
    Indel_rr$VAF <-  as.numeric(geno(Indels_Sample)$VAF)
    Indel_rr$FILTER[Indel_rr$VAF < 0.25] <- "FAIL_VAF"
    Indel_rr <- get_indel_context(Indel_rr, ref_genome)
    
    indels_raw[[Sample]] <- Indel_rr
  }
}

### Read the blacklist
exclusion_vcf <- readVcf(paste(Input_dir, "Indels/excludelist_v1.vcf", sep = ""), "hg38")
exclusion_grl <- rowRanges(exclusion_vcf)
exclusion_grl <- exclusion_grl[seqnames(exclusion_grl) %in% c(1:22, "X", "Y"),]
seqlevels(exclusion_grl) <- seqlevels(exclusion_grl)[1:24] # The decoy chromosomes give errors later on
seqlevels(exclusion_grl) <- paste("chr", seqlevels(exclusion_grl), sep = "")
exclusion_grl <- get_mut_type(exclusion_grl, type = "indel") # there appears to be one SNVs according to MutationalPatterns

### Characterize the indels in the blacklist
### !!! This takes a long time
#blacklist_grl <- get_indel_context(blacklist_grl, ref_genome)
exclusion_context <- GRanges()
for(chr in seqlevels(exclusion_grl)){
  print(chr)
  exclusion_chr <- exclusion_grl[which(seqnames(exclusion_grl) == chr)]
  exclusion_chr_context <- get_indel_context(exclusion_chr, ref_genome)
  exclusion_context <- c(exclusion_context, exclusion_chr_context)
}

exclusion_counts <- count_indel_contexts(exclusion_context)
colnames(exclusion_counts) <- "PTATO\nIndel\nExclusion\nList"

# Plot how often each indel in the list occurs
RecurrencyList_Info <- read.table(paste(Input_dir, "Indel_RecurrencyList_Info.txt", sep = ""), sep = ";")

RecurrencyList_Info$V4 <- gsub(RecurrencyList_Info$V4, pattern = "IF=", replacement = "")
RecurrencyList_Info$V3 <- gsub(RecurrencyList_Info$V3, pattern = "IC=", replacement = "")
RecurrencyList_Info$V2 <- gsub(RecurrencyList_Info$V2, pattern = "SF=", replacement = "")

ggplot(RecurrencyList_Info, aes(x = as.numeric(V3))) + geom_histogram(bins = 21, fill = "dodgerblue4", col = "white", linewidth = 0.1) + labs(y = "Indels in exclusion list (#)", x = "Individuals (#)") +
  theme_classic(base_size = 6) + scale_x_continuous(breaks = seq(0,25,2)) + scale_y_continuous(breaks = seq(0,2e6, 2.5e5), expand = c(0,10000))
ggsave(paste(Output_dir_figS2, "FigureS2b_Indel_Histogram.pdf", sep = ""), width = 7, height = 4, units = "cm")

#ggplot(RecurrencyList_Info, aes(x = as.numeric(V3))) + geom_bar(stat = "count")
ggplot(RecurrencyList_Info[1:10000,], aes(x = as.numeric(V4), y = as.numeric(V2))) + geom_point()

# # Determine how well indel filtering works if indels are found in 3 or more individuals (instead of 2+)
# exclusion3ind <- exclusion_vcf[which(RecurrencyList_Info$V3 > 2),]
# exclusion3ind_grl <- rowRanges(exclusion3ind)
# exclusion3ind_grl <- exclusion3ind_grl[seqnames(exclusion3ind_grl) %in% c(1:22, "X", "Y"),]
# seqlevels(exclusion3ind_grl) <- seqlevels(exclusion3ind_grl)[1:24] # The decoy chromosomes give errors later on
# seqlevels(exclusion3ind_grl) <- paste("chr", seqlevels(exclusion3ind_grl), sep = "")
# exclusion3ind_grl <- get_mut_type(exclusion3ind_grl, type = "indel") # there appears to be one SNVs according to MutationalPatterns

# Filter the indels for blacklist and homopolymers
indel_grl_filter <- filter_indels(indel_grl = indels_raw, indel_recurrency_grl = exclusion_grl, remove_homopolymer = TRUE)
#indel_grl_filter3 <- filter_indels(indel_grl = indels_raw, indel_recurrency_grl = exclusion3ind_grl, remove_homopolymer = TRUE)

indel_grl <- lapply(indel_grl_filter, function(x) x[which(x$FILTER =="PASS"),])
#indel_grl3 <- lapply(indel_grl_filter3, function(x) x[which(x$FILTER =="PASS"),])

indels_grl_raw <- lapply(indels_raw, function(x) x[which(x$FILTER =="PASS"),])

## 1g: Indel spectra before and after filtering
indel_counts_raw <- count_indel_contexts(indels_grl_raw)
indel_counts <- count_indel_contexts(indel_grl)

## Adapted from MutationalPatterns
# Reformat to long format
indel_counts_m <- indel_counts %>% as.data.frame() %>% tibble::rownames_to_column("muttype_total") %>% 
  tidyr::separate(muttype_total, c("muttype", "muttype_sub"), 
                  sep = "_(?=[0-9])") %>% dplyr::mutate(muttype = factor(muttype, 
                                                                         levels = unique(muttype)))
indel_counts_m <- indel_counts_m %>% dplyr::select(-muttype_sub) %>% 
  dplyr::group_by(muttype) %>% dplyr::summarise_all(list(~sum(.))) %>% 
  tidyr::gather(key = "Sample", value = "count", -.data$muttype) %>% 
  dplyr::mutate(Sample = factor(Sample, levels = unique(Sample)))
indel_counts_m$Filter <- "After PTATO"

indel_counts_raw_m <- indel_counts_raw %>% as.data.frame() %>% tibble::rownames_to_column("muttype_total") %>% 
  tidyr::separate(muttype_total, c("muttype", "muttype_sub"), 
                  sep = "_(?=[0-9])") %>% dplyr::mutate(muttype = factor(muttype, 
                                                                         levels = unique(muttype)))
indel_counts_raw_m <- indel_counts_raw_m %>% dplyr::select(-muttype_sub) %>% 
  dplyr::group_by(muttype) %>% dplyr::summarise_all(list(~sum(.))) %>% 
  tidyr::gather(key = "Sample", value = "count", -.data$muttype) %>% 
  dplyr::mutate(Sample = factor(Sample, levels = unique(Sample)))
indel_counts_raw_m$Filter <- "Before PTATO"

# Merge the raw and the filtered indel_counts
indel_counts_merged <- rbind(indel_counts_m, indel_counts_raw_m)

# Add the names of the individuals and the (shortened) sample labels
indel_counts_merged <- merge(indel_counts_merged, Metadata[,c("Sample", "Individual", "Label")], all.x = T)
#indel_counts_merged$label <- indel_counts_merged$count
#indel_counts_merged$label[indel_counts_merged$count < 50 & indel_counts_merged$Filter == "Raw"] <- ""
#indel_counts_merged$label[indel_counts_merged$count < 5 & indel_counts_merged$Filter == "Filtered"] <- ""

indel_counts_merged$Filter <- factor(indel_counts_merged$Filter, levels = c("Before PTATO", "After PTATO"))
# Shorten the legend labels
levels(indel_counts_merged$muttype) <- gsub(pattern = "microhomology", replacement = "MH", x =levels(indel_counts_merged$muttype)  )
levels(indel_counts_merged$muttype) <- gsub(pattern = "deletion", replacement = "del", x =levels(indel_counts_merged$muttype)  )
levels(indel_counts_merged$muttype) <- gsub(pattern = "insertion", replacement = "ins", x =levels(indel_counts_merged$muttype)  )

indel_counts_merged$Label <- factor(indel_counts_merged$Label, levels =  c("AMLBULK","AML-PTAD2B1", "AML-PTAD2D1" , "HSC-PTAP1B8","HSC-PTAP1D9" ,"HSC-PTAP1E9" ,"HSC-PTAP1G9" ,"HSP1K6",  "HSP3B22"  , "HSP3D22",  "MPP-PTAP6F4" ,"MPP-PTAP6F5" ,"MPP-PTAP6F6"))

indel_colors <- c("#FDBE6F", "#FF8001", "#B0DD8B", "#36A12E", "#FDCAB5", 
                  "#FC8A6A", "#F14432", "#BC141A", "#D0E1F2", "#94C4DF", 
                  "#4A98C9", "#1764AB", "#E2E2EF", "#B6B6D8", "#8683BD", 
                  "#61409B")

ggplot(indel_counts_merged, aes(x = Label, y = count, 
                                #label = as.numeric(label), 
                                fill = muttype)) + 
  geom_bar(stat = "identity",  col = "black", width = 0.75, linewidth = 0.2) + 
  
  #geom_text( size = 1.8, position = position_stack(vjust = 0.5)) +
  geom_text(aes(label = round(after_stat(y), 0), group = Label), 
            stat = 'summary', fun = sum, vjust = -0.3, size = 1.8) +
  facet_grid(Filter ~ Individual, scales = "free", space = "free_x", switch = "y") +
  labs(x = "Sample", y = "Somatic Indels (#)", fill = "Indel type") + 
  scale_fill_manual(values = indel_colors) + 
  scale_y_continuous(expand = expansion(mult = c(0,0.1))) +   
  theme_bw(base_size = 6) + 
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1),
        legend.key.height = unit(0.3, 'cm'),
        legend.key.width = unit(0.3, 'cm'),
        panel.spacing.x  = unit(0, "lines"), strip.placement = "outside",
        legend.position="bottom")

ggsave(paste(Output_dir_fig1, "Fig1h_Indels_PTATO.pdf", sep = ""), width = 7.5, height = 7.5, units = "cm")

# plot the full indel spectra for two IBFM35 samples as an example
indel_counts_IBFM35 <- cbind(indel_counts_raw[,c("IBFM35-DX2BM-AMLBULK","IBFM35-DX2BM-HSCPTAP1E9")],indel_counts[,"IBFM35-DX2BM-HSCPTAP1E9"])
colnames(indel_counts_IBFM35) <- c("IBFM35\nAMLBULK", "HSC-PTAP1E9\nBefore PTATO", "HSC-PTAP1E9\nAfter PTATO")

pdf(file = paste(Output_dir_fig1, "Fig1g_Indels_IBFM35.pdf", sep = ""), width = 4, height = 2.75, onefile = F)
plot_indel_contexts2(counts = indel_counts_IBFM35, condensed = T)
dev.off()

### Blacklist content

pdf(file = paste(Output_dir_FigS2, "FigureS2a_Indel_context_ExclusionList.pdf", sep = ""), width = 5, height = 2, onefile = F)
plot_indel_contexts2(counts = blacklist_counts, condensed = T)
dev.off()

### Determine the reason why indels are removed from samples 
indels_flagged <- filter_indels(indel_grl = indels_grl_raw, indel_recurrency_grl = exclusion_grl, remove_homopolymer = TRUE, mode = "flag")
indels_filtered <- filter_indels(indel_grl = indels_grl_raw, indel_recurrency_grl = exclusion_grl, remove_homopolymer = TRUE)

indel_counts <- data.frame(Sample = names(indels_filtered), 
                           PASS = lengths(indels_filtered), 
                           RECURRENT = lengths(lapply(indels_flagged, function(x) x[which(x$FILTER =="RECURRENT"),]) ),
                           HOMOPOLYMER = lengths(lapply(indels_flagged, function(x) x[which(x$FILTER =="HOMOPOLYMER"),]) ),
                           RECURRENT_HOMOPOLYMER = lengths(lapply(indels_flagged, function(x) x[which(x$FILTER =="RECURRENT;HOMOPOLYMER"),]) ))

indel_counts$Individual <- gsub(pattern = "-.*", replacement = "", indel_counts$Sample)
indel_counts$Label <- gsub(pattern = ".*-", replacement = "", indel_counts$Sample)

indel_counts_m <- melt(indel_counts[,c("Label", "Individual", "PASS", "RECURRENT", "HOMOPOLYMER", "RECURRENT_HOMOPOLYMER")])

indel_counts_m$variable <- factor(indel_counts_m$variable, levels = c("RECURRENT", "RECURRENT_HOMOPOLYMER", "HOMOPOLYMER", "PASS"))
# Bulk samples are indicated with an asterisk
indel_counts_m$Label[!grepl("PTA", indel_counts_m$Label)] <- paste(indel_counts_m$Label[!grepl("PTA", indel_counts_m$Label)], "*", sep = "")
ggplot(indel_counts_m, aes(x = Label, y = value, fill = variable)) + 
  geom_bar(stat = "identity", width = 0.8) +
  facet_grid(~Individual, space = "free", scales = "free") +
  scale_fill_manual(values = c("#D55E00","#a6761d", "#E6AB02","#009E73")) +
  scale_y_continuous(expand = c(0,0), limit= c(0,1250), breaks = seq(0,1250, 250)) +
  theme_classic(base_size = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.25, 'cm')) +
  labs(x = "Sample", y = "Indels (VAF>0.25)", fill = "Filter")

ggsave(paste(Output_dir_FigS2, "FigureS2c_Indel_Filter.pdf", sep = ""), width = 11, height = 5, units = "cm")

# plot the indel profile of the exclusion list
recurrency_indel_counts <- count_indel_contexts(indel_recurrency_grl)
colnames(recurrency_indel_counts) <- "Recurrent\nIndel list"
pdf(file = paste(Output_dir_ExtData2, "Figure2Sb_Indel_ExclusionList.pdf", sep = ""), width = 5, height = 2, onefile = F)
plot_indel_contexts2(counts = recurrency_indel_counts, condensed = T, y_label = " (#)")
dev.off()


### Figure S2d: Cosine similarities indel spectra bulk/PTA before and after PTATO filtering
bulk_samples <- Metadata$Sample[Metadata$Type == "Bulk"]
PTA_samples <- Metadata$Sample[Metadata$Type != "Bulk" & Metadata$Individual != "PMCCB15"]

makeIndelMatrix <- function(Samples, indel_counts, filter){
  indel_matrix <- data.frame()
  for(Sample in Samples){
    indels_sample <- indel_counts[indel_counts$Sample == Sample & indel_counts$Filter == filter,]
    print(Sample)
    if(nrow(indel_matrix) == 0){
      indel_matrix <- indels_sample[,c("muttype", "count")]
      colnames(indel_matrix)[ncol(indel_matrix)] <- Sample
      
    } else {
      indel_matrix <- merge(indel_matrix, indels_sample[,c("muttype", "count")], by = "muttype")
      colnames(indel_matrix)[ncol(indel_matrix)] <- Sample
    }
  }
  
  indel_matrix2 <- as.matrix(indel_matrix[,-1])
  row.names(indel_matrix2) <- indel_matrix[,1]
  return(indel_matrix2)
}

Indels_bulk_before <- makeIndelMatrix(Samples = bulk_samples, indel_counts = indel_counts_merged, filter = "Before PTATO")
Indels_bulk_after <- makeIndelMatrix(Samples = bulk_samples, indel_counts = indel_counts_merged, filter = "After PTATO")

Indels_PTA_before <- makeIndelMatrix(Samples = PTA_samples, indel_counts = indel_counts_merged, filter = "Before PTATO")
Indels_PTA_after  <- makeIndelMatrix(Samples = PTA_samples, indel_counts = indel_counts_merged, filter = "After PTATO")

Indels_merged <- data.frame(Bulk_before = Indels_bulk_before, PTA_before = Indels_PTA_before, Bulk_after = Indels_bulk_after, PTA_after = Indels_PTA_after)
cos_sim_m = cos_sim_matrix(Indels_merged, Indels_merged)

#plot_cosine_heatmap(cos_sim_m, cluster_rows = F, cluster_cols = F, plot_values = T)

Indels_summed <- data.frame(Bulk_before = rowSums(Indels_bulk_before), PTA_before = rowSums(Indels_PTA_before), Bulk_after = rowSums(Indels_bulk_after), PTA_after = rowSums(Indels_PTA_after))
Indels_summed <- as.matrix(Indels_summed)

# Indels_after_merged <- data.frame(Bulk = rowSums(Indels_bulk_after), PTA = rowSums(Indels_PTA_after))
# Indels_after_merged <- as.matrix(Indels_after_merged)

cos_sim_matrix = cos_sim_matrix(Indels_summed, Indels_summed)
cos_sim_matrix_m <- melt(cos_sim_matrix)
cos_sim_matrix_m$Var1 <- gsub(cos_sim_matrix_m$Var1, pattern = "_", replacement = "\n")
cos_sim_matrix_m$Var2 <- gsub(cos_sim_matrix_m$Var2, pattern = "_", replacement = "\n")

cos_sim_matrix_m$Var2 <- factor(cos_sim_matrix_m$Var2, levels = c( "PTA\nafter","Bulk\nafter" ,"PTA\nbefore", "Bulk\nbefore"))

cos_sim_matrix_m$Var1 <- factor(cos_sim_matrix_m$Var1, levels = c("Bulk\nbefore","PTA\nbefore", "Bulk\nafter" , "PTA\nafter"))
cos_sim_matrix_m$Var2
# cos_sim_matrix_m$facet_a <- ifelse(grepl(cos_sim_matrix_m$Var1, pattern = "before") == T, "Before", "After")
# cos_sim_matrix_m$facet_b <- ifelse(grepl(cos_sim_matrix_m$Var2, pattern = "before") == T, "Before", "After")
# cos_sim_matrix_m$facets <- factor(paste(cos_sim_matrix_m$facet_a, cos_sim_matrix_m$facet_b, sep = "_"))

ggplot(cos_sim_matrix_m, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile(col = "white", linewidth = 0.5) +
  geom_text(aes(label = round(value, 2)), size = 2) + 
  scale_x_discrete(position = "top")  +
  scale_fill_gradient(
    low = "cornsilk",
    high = "#1F78B4",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill") +
  labs(fill = "Cosine\nSimilarity", x = "Main indel profile", y = "Main indel profile") +
  theme_classic(base_size = 6) +
  theme(        strip.background = element_blank(),
                strip.text.x = element_blank(),
                legend.key.height= unit(0.3, 'cm'),
                legend.key.width= unit(0.3, 'cm'), panel.border=element_blank(), axis.line=element_line())

ggsave(paste(Output_dir_FigS2, "FigureS2d_Indel_InputOutput_CosineSim.pdf", sep =""), width = 7, height = 5.5, units = "cm")




