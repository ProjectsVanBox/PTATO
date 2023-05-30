# Supplemental Figures S1
library(ggplot2)
library(VariantAnnotation)
library(MutationalPatterns)
library(reshape2)
library(pheatmap)
library(ggh4x)
library(randomForest)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)

# Path to the Somatic_Data directory (eg "/PTATO/Somatic_Data/" )
Input_dir <-  "/path/to/MendeleyData/directory/"
# Path to the directory containing the GeneralFunctions.R and manual_training_performance_functions.R scripts (eg "/PTATO/Scripts/" )
Scripts_dir <- "/path/to/GithubScripts/directory/"

Output_dir_FigS1 <- "/path/to/output/directory/"
if(dir.exists(Output_dir_FigS1) == F){dir.create(Output_dir_FigS1)}

source(paste(Scripts_dir, "manual_training_performance_functions.R", sep = ""))
source(paste(Scripts_dir, "GeneralFunctions.R", sep = ""))

Metadata_file <- paste(Input_dir, "Table_S1.txt", sep = "")
Metadata <- read.delim(Metadata_file, header = T)

ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

### Figure S1b: Comparison of performance of RF models trained on different ratios of true and false positives
cutoff_ratio <- matrix(ncol=13)[-1,]

for ( x in c("0.2_0.8","0.3_0.7","0.4_0.6","0.5_0.5","0.6_0.4","0.7_0.3","0.8_0.2","0.9_0.1")) {
  print(x)
  res <- readRDS(paste(Input_dir, "/Training/Manual/",x,"/res.rds",sep=""))
  output_forest <- readRDS(paste(Input_dir, "/Training/Manual/",x,"/randomforest_v3.0.1.rds",sep=""))
  
  res$pta_prob <- output_forest$all$votes[,"PTA"]
  
  cutoff_accuracies = lapply(seq(0,1,0.01), determine_cutoff_accuracy, res) %>%
    bind_rows() %>%
    dplyr::mutate(precision = ifelse(is.nan(precision), 1 ,precision),
                  F1 = ifelse(is.nan(F1), 0, F1))
  cutoff_accuracies$ratio = x
  
  cutoff_ratio <- rbind(cutoff_ratio, cutoff_accuracies)
}  

cutoff_ratio$avg_class_error <- (cutoff_ratio$noPTA_class_error+cutoff_ratio$PTA_class_error)/2

cutoff_ratio$ratio <- gsub(cutoff_ratio$ratio, pattern = "_", replacement = " : ")

ggplot(data=cutoff_ratio, aes(x=cutoff,y=1-avg_class_error,col=ratio)) +
  geom_line() + theme_classic(base_size = 6) +
  theme(legend.spacing.y = unit(0.2, 'cm'), legend.key.height = unit(0.2, "cm")) +
  labs(y = "Balanced Accuracy", x = "PTA probability cutoff", col = "Ratio in training set\n(True pos : Artefacts)")
ggsave(paste(Output_dir_FigS1, "FigS1A_Training_TP_FP_Ratios.pdf", sep = ""), width = 8, height = 4, units = "cm")

# Figure S1c: changes in accuracy when iteratively removing features
feature_remove_dir <- paste(Input_dir, "Training/Manual/LooFeatures/", sep = "")
meandecreasegini <- read.table(paste(Input_dir, "Training/Manual/LooFeatures/randomforest_v3.0.2_importance.txt", sep = ""), header = T)
meandecreasegini <- meandecreasegini[order(meandecreasegini[,1], decreasing = F),, drop = F]
feature_overview <- data.frame()
for(i in 1:25){
  print(i)
  feature_data_file <- list.files(path = paste(feature_remove_dir, i, "/", sep = ""), pattern = "_context.txt", full.names = T, recursive = T)
  feature_data <- read.table(feature_data_file, skip = 25, nrows = 1)
  
  feature_overview <- rbind(feature_overview, data.frame(Removed = i, BalancedAccuracy = feature_data[,3]))
}

feature_overview$Features <- 26-feature_overview$Removed
feature_overview$Feature <- row.names(meandecreasegini)[1:25]
feature_overview$Label <- paste(feature_overview$Removed, " (-",feature_overview$Feature, ")", sep = "")
#ggplot(feature_overview, aes( x= Removed, y = BalancedAccuracy)) + geom_line()
ggplot(feature_overview, aes( x= Removed, y = BalancedAccuracy)) + 
  geom_line(col = "#7570b3") + 
  coord_cartesian(ylim = c(0.5, 0.75)) + 
  theme_classic(base_size = 6) + labs(x = "Features excluded for training (#)", y = "Balanced accuracy\nof RF model") + 
  scale_x_continuous(breaks = seq(1,25,1),labels = feature_overview$Label ) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste(Output_dir_FigS1, "FigureS1B_Feature_accuracy_drop.pdf", sep = ""), width = 7, height = 4, units = "cm")

# Read the PTATO SNV VCFs
PTATO_Data <- list() 
snv_grl_raw <- list()
for(Individual in c("PMCCB15", "IBFM35")){
  print(Individual)
  
  PTATO_snv_dir <- paste(Input_dir, Individual, "/PTATO/snvs/", sep = "")
  print(PTATO_snv_dir)
  PTATO_indel_dir <- paste(Input_dir, Individual, "/PTATO/indels/", sep = "")
  
  
  for(Sample in Metadata$Sample[Metadata$Individual == Individual & Metadata$Germline == FALSE]){
    print(Sample)
    if(Sample != "PMCCB15-CBWTARAC-HSP1K6"){
      PTATO_Data[[Individual]][[Sample]][["vcf_unfiltered"]] <- list.files(PTATO_snv_dir, recursive = T, pattern = ".ptato.vcf.gz$", full.names = T)[grep(Sample, x =  list.files(PTATO_snv_dir, recursive = T, pattern = ".ptato.vcf.gz$", full.names = T))]
      PTATO_Data[[Individual]][[Sample]][["vcf_filtered"]] <- list.files(PTATO_snv_dir, recursive = T, pattern = ".ptato.filtered.vcf.gz$", full.names = T)[grep(Sample, x =  list.files(PTATO_snv_dir, recursive = T, pattern = ".ptato.filtered.vcf.gz$", full.names = T))]
      
      snv_grl_raw[[Sample]] <-  readPTATOvcf(vcf_unfiltered = PTATO_Data[[Individual]][[Sample]][["vcf_unfiltered"]], 
                                             vcf_filtered = PTATO_Data[[Individual]][[Sample]][["vcf_filtered"]], 
                                             VAF_threshold = 0.2)

    }
  }
}

# Color scheme
Sample_colors <- c(brewer.pal(n = 7, name = "Dark2"))
names(Sample_colors) <- c(Metadata$Label[Metadata$Individual %in% c("PMCCB15", "IBFM35") & Metadata$Training == "Yes"])

### Figure S1c: Plot PTAprobs distribution
# !!! ADD CLONE TO PLOTS
for(Individual in c("PMCCB15", "IBFM35")){
  print(Individual)
  Probs <- data.frame()
  for(Sample in Metadata$Sample[Metadata$Individual == Individual & Metadata$Germline == FALSE & Metadata$Training == "Yes"]){
    print(Sample)
    Probs_sample <- data.frame(Sample = Sample, 
                               ID = Metadata$Label[Metadata$Sample == Sample],
                               Type = Metadata$Type[Metadata$Sample == Sample],
                               PTAprob = snv_grl_raw[[Sample]]$PTAprob, 
                               PTAprobCutoff = snv_grl_raw[[Sample]]$PTAprobCutoff, 
                               TotalVariants = length(snv_grl_raw[[Sample]]))
    Probs <- rbind(Probs, Probs_sample)
  }
  Cutoffs <- Probs[!duplicated(Probs$ID), c("ID","PTAprobCutoff", "Type")]
  names(Cutoffs)[2] <- "PTAprob"
  ggplot(Probs, aes(x = PTAprob, color = ID)) + 
    stat_density( aes(linetype = Type),
                 geom="line",position="identity", size = 0.3) + 
    #geom_density(size = 0.3, aes(linetype = Type)) + 
    geom_vline(data = Cutoffs, aes(xintercept = as.numeric(PTAprob), color = ID), alpha = 0.7, linetype = 2, size = 0.3, show.legend=FALSE ) +
    ggtitle(Individual) +
    scale_y_continuous(expand = c(0, 0))  + 
    scale_x_continuous(expand = c(0, 0)) +
    scale_color_manual(values = Sample_colors, breaks = Probs$ID) +
    scale_linetype_manual(values=c("Bulk" = 5, "PTA" = 1), guide = "none")+
    theme_classic(base_size = 6) + 
    theme(panel.grid.major = element_blank(),plot.title = element_text(hjust = 0.5),
          legend.spacing.y = unit(0.2, 'cm'),
          legend.key.height = unit(0.3, 'cm'),
          legend.key.width = unit(0.5, 'cm'), legend.position = "bottom") +
    labs(col = "Sample")
  ggsave(paste(Output_dir_FigS1, Individual, "_PTAprobDensity.pdf", sep = ""), width = 4, height = 4, units = "cm")
  
  ggplot(Probs, aes(x = PTAprob, col = ID)) +
    geom_step(aes(len=TotalVariants,y=..y.. * len), stat="ecdf", size = 0.3) + 
    geom_vline(aes(xintercept = as.numeric(PTAprobCutoff), col = ID), linetype = 2, alpha = 0.7, size = 0.3, show.legend=FALSE ) +
    ggtitle(Individual) +
    scale_y_continuous(expand = c(0, 0))  + 
    scale_x_continuous(expand = c(0, 0)) +
    scale_color_manual(values = Sample_colors, breaks = Probs$ID) +
    scale_linetype_manual(values=c("Bulk" = 5, "PTA" = 1), guide = "none")+
    theme_classic(base_size = 6) + 
    theme(panel.grid.major = element_blank(),plot.title = element_text(hjust = 0.5),
          legend.spacing.y = unit(0.2, 'cm'),
          legend.key.height = unit(0.3, 'cm'),
          legend.key.width = unit(0.5, 'cm'), legend.position = "bottom") +
    labs(col = "Sample", y = "Base substitutions (#)")
  ggsave(paste(Output_dir_FigS1, Individual, "_PTAprobCumulative.pdf", sep = ""), width = 4, height = 4, units = "cm")
}

## Figure S1d: Precision-recall curves
precision_merged <- data.frame()

for(Individual in c("PMCCB15", "IBFM35")){
  print(Individual)
  
  PTATO_snv_dir <- paste(Input_dir, Individual, "/PTATO/snvs/", sep = "")
  print(PTATO_snv_dir)
  ptatotable_files <- list.files(path = PTATO_snv_dir, pattern = ".ptatotable.txt", recursive = T, full.names = T)
  
  for(Sample in Metadata$Sample[Metadata$Individual == Individual & Metadata$Germline == FALSE & Metadata$Training == "Yes"]){
    print(Sample)
    ptatotable_file <- ptatotable_files[grep(Sample, ptatotable_files)]
    print(ptatotable_file)
    
    ptatotable <- read.table(ptatotable_file, header = T)
    
    ptatotable$prec_recall <- abs(ptatotable$precision-ptatotable$recall)
    
    ptatotable2 <- ptatotable[ptatotable$precision > 0 & ptatotable$recall > 0,]
    
    lowest <- ptatotable2[which(ptatotable2$prec_recall == min(ptatotable2$prec_recall, na.rm = T)),]
    
    precision <- melt(ptatotable[,c("PTAprob_cutoff", "precision", "recall")], id.vars = "PTAprob_cutoff")

    levels(precision$variable) <- c("Precision", "Recall")
    precision$Sample <- Sample
    precision$Cutoff <- mean(lowest$PTAprob_cutoff)
    precision_merged <- rbind(precision_merged, precision)
    
    
    ggplot(precision, aes(x = PTAprob_cutoff, y = value, col = variable)) + 
      geom_line(size = 0.4) + 
      geom_vline(aes(xintercept = as.numeric(Cutoff)), linetype= 5, size = 0.3, color = "#999999") +
      ggtitle(Sample) +
      theme_classic(base_size = 6) + 
      scale_y_continuous(expand = c(0, 0.01)) + 
      scale_x_continuous(expand = c(0, 0)) +
      scale_color_manual(values = c("#0072B2", "#D55E00")) +
      theme(legend.position = c(0.75, 0.1),
            panel.border = element_rect(colour = "black", fill=NA, size = 0.5),
            legend.title=element_blank(),
            legend.background = element_rect(fill=alpha("white", 0.9)),
            legend.spacing.y = unit(0, 'cm'),
            legend.key.height = unit(0.1, 'cm'),
            legend.key.width = unit(0.5, 'cm'),
            plot.title = element_text(hjust = 0.5, size = 6)) +
      guides(color = guide_legend(byrow = TRUE)) + 
      labs(y = "Value", col = "Variable", x = "PTAprob Linked-read Cutoff")
    ggsave(paste(Output_dir_FigS1, Sample, "_PrecisionRecall.pdf", sep = ""), width = 4, height = 3.5, units = "cm")
    
  }
}

## Figure S1e: Precision-recall values of all samples in training set
mean_PrecisionRecall <- data.frame()
for(Sample in unique(precision_merged$Sample)){
  
  precision_sample <- precision_merged[precision_merged$Sample == Sample,]
  mean_PrecisionRecall_sample <- data.frame(Sample = Sample, 
                                            Label = Metadata$Label[Metadata$Sample == Sample], 
                                            Individual =  Metadata$Individual[Metadata$Sample == Sample],
                                            Type =  Metadata$Type[Metadata$Sample == Sample],
                                            PrecisionRecall = mean(precision_sample$value[precision_sample$PTAprob_cutoff == round(unique(precision_sample$Cutoff), 2)]))
  mean_PrecisionRecall <- rbind(mean_PrecisionRecall, mean_PrecisionRecall_sample)
}

ggplot(mean_PrecisionRecall, aes(x = Label, y = PrecisionRecall, fill = Type)) + 
  geom_bar(stat = "identity", show.legend = F, width = 0.8) +
  facet_grid(~Individual, scales = "free", space = "free") +
  theme_grey(base_size = 6) +
  labs(y = "Precision/Recall\nat linked-read Cutoff", x = "Sample") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + 
  scale_fill_manual(values = c("darkgray", "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        #strip.placement =  "outside",
        strip.background = element_rect(color = "black" , size = 0.3),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.3))
ggsave(paste(Output_dir_FigS1, "FigureS1F_PrecisionRecall_Overview.pdf", sep = ""), width = 5, height = 3.5, units = "cm")

### 

mut_mat_sample <- calc_cosim_mutmat(gr = snv_grl_raw[["PMCCB15-CBMPP-PTAP6F4"]])
cos_sim_matrix <- cos_sim_matrix(mut_mat_sample[["PASS"]], mut_mat_sample[["PASS"]])

### CLustering (code from MutationalPatterns)
hc.sample <- hclust(dist(cos_sim_matrix[7:nrow(cos_sim_matrix),
                                        7:ncol(cos_sim_matrix)]), method = "complete")

# This function can be used to get the samples in two clusters
sub_grp <- cutree(hc.sample, k = 2)
PTAprob_cutoff <- max(names(sub_grp[sub_grp == 1]))
sub_grp_df <- data.frame(Cluster = sub_grp)
sub_grp_df$Cluster[sub_grp_df$Cluster == 1] <- "Low"
sub_grp_df$Cluster[sub_grp_df$Cluster == 2] <- "High"
annot_colors=list(Cluster=c(High="#7570b3",Low="#1b9e77"))
pdf(file = paste(Output_dir_FigS1,"FigureS1G_CosineHeatmap.pdf", sep = ""), width = 4, height = 3, pointsize = 5)

pheatmap(cos_sim_matrix[7:nrow(cos_sim_matrix),
                        7:ncol(cos_sim_matrix)],
         annotation_row = sub_grp_df, annotation_col = sub_grp_df,
         cutree_rows = 2,
         cutree_cols = 2,
         color = colorRampPalette(brewer.pal(n = 7, name =
                                                   "YlOrRd"))(100),
         annotation_colors=annot_colors, border_color = NA, fontsize = 5)
dev.off()

# Plot the 96-trinucleotide profiles of the variants for two different cutoffs (as an example)
snv_03 <- snv_grl_raw[["PMCCB15-CBMPP-PTAP6F4"]][which(snv_grl_raw[["PMCCB15-CBMPP-PTAP6F4"]]$PTAprob <= 0.325),]
snv_07 <- snv_grl_raw[["PMCCB15-CBMPP-PTAP6F4"]][which(snv_grl_raw[["PMCCB15-CBMPP-PTAP6F4"]]$PTAprob <= 0.7),]

mut_mat_comparison <- cbind(mut_matrix(snv_03,ref_genome = ref_genome), mut_matrix(snv_07,ref_genome = ref_genome) )
colnames(mut_mat_comparison) <- c("Cutoff = 0.325", "Cutoff = 0.7")
pdf(paste(Output_dir_FigS1, "FigureS1G_96profiles.pdf", sep = ""), width = 3, height = 2, onefile = F)

plot_96_profile_small(mut_mat_comparison, condensed = T, ymax = 0.1)
dev.off()

# Print the cosine similarity between the variants above and below the cutoffs
cos_sim_matrix(mut_mat_comparison,mut_mat_comparison)
