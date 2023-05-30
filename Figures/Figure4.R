# Figures 4 and 5I and S9C

library(MutationalPatterns)
library(ggplot2)
library(reshape2)
library(VariantAnnotation)
library(RColorBrewer)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
library(nlme)
library(ggeffects)
library(lme4)
# For facet_nested:
library(ggh4x)
# For statistics:
library(ggpubr)
library(rstatix)

Input_dir <- "/path/to/MendeleyData/directory/"
Output_dir <-  "/path/to/output/directory/"
Output_dir_FigureS9 <- "/path/to/output/directory/"
Scripts_dir <- "/path/to/GithubScripts/directory/"

if(dir.exists(Output_dir) == F){
  print(paste("# Creating output dir:",Output_dir, sep =""))
  dir.create(Output_dir, recursive = T)
}

source(paste(Scripts_dir, "GeneralFunctions.R", sep = ""))

Metadata <- read.delim(paste(Input_dir, "Table_S1.txt", sep = ""))
Metadata_FA <- Metadata[Metadata$Phenotype == "FanconiAnemia" | Metadata$Individual == "IBFM35",]
Metadata_FA <- Metadata_FA[Metadata_FA$Germline == FALSE,]
Metadata_Ageline <- Metadata[Metadata$Phenotype == "Healthy" & Metadata$Training == "No" & Metadata$Type == "Clone" & Metadata$Individual != "PMCCB15" & Metadata$Individual != "STE0072",]
Metadata_Fig4 <- rbind(Metadata_FA,Metadata_Ageline)

### Read the PTATO vcfs
Indels_raw <- list()
MinimalVAF <- 0.25 # all variants will have a VAF higher than MinimalVAF
for(Sample in Metadata_FA$Sample){
  print(Sample)
  # For PTA samples we use the PTATO-filtered VCFs, for bulk samples we use the SMuRF-filtered VCFs
  if(Metadata_FA$Type[Metadata_FA$Sample == Sample] == "PTA"){
    # PTATO Indels
    indels_vcf_file <- list.files(paste(Input_dir, Metadata_FA$Individual[Metadata_FA$Sample == Sample], "/PTATO/", sep = ""),
                                  pattern = paste(Sample, ".indels.ptato.callable.vcf.gz$", sep = ""), full.names = T, recursive = T)
    print(indels_vcf_file)
    if(file.exists(indels_vcf_file) == T){
      Indels_sample <- readVcf(indels_vcf_file, genome = "hg38")
      Indels_sample <- Indels_sample[as.vector(seqnames(Indels_sample)) %in% c(1:22, "X", "Y"),]
      Indels_raw[[Sample]] <- VCF_to_GR(Indels_sample)
      Indels_raw[[Sample]]$VAF <- geno(Indels_sample)$VAF
      Indels_raw[[Sample]] <- Indels_raw[[Sample]][which(Indels_raw[[Sample]]$VAF > MinimalVAF)]
    } else {
      stop("!!! VCF not found")
    }

  } else if(Metadata_FA$Type[Metadata_FA$Sample == Sample]%in% c("Clone","Bulk")){
    # SMuRF VCFs contain both somatic SBSs and indels
    vcf_file <- list.files(paste(Input_dir, Metadata_FA$Individual[Metadata_FA$Sample == Sample], "/PTATO/", sep = ""),
                           pattern = paste(Sample, ".SMuRF.filtered.callable.vcf.gz$", sep = ""), full.names = T, recursive = T)
    print(vcf_file)
    vcf <- readVcf(vcf_file, genome = "hg38")
    vcf <- vcf[as.vector(seqnames(vcf)) %in% c(1:22, "X", "Y"),]
    Indels_sample <- vcf[isIndel(vcf),]
    Indels_raw[[Sample]] <- VCF_to_GR(Indels_sample)
    Indels_raw[[Sample]]$VAF <- geno(Indels_sample)$VAF
    Indels_raw[[Sample]] <- Indels_raw[[Sample]][which(Indels_raw[[Sample]]$VAF > MinimalVAF)]
    
  }
}

# Only autosomal
Indels_raw_autosomal <-  lapply(Indels_raw, function(x) x[which(as.vector(seqnames(x)) %in% paste("chr", c(1:22), sep = "")),])
Indels_raw_autosomal <- get_mut_type(Indels_raw_autosomal, type = "indel")
Indels_raw_autosomal <- get_indel_context(Indels_raw_autosomal, ref_genome)



### Read the indel filtering list
indel_exclusion_vcf <- readVcf(paste(Input_dir,"Indels/excludelist_v1.vcf", sep = ""),"hg38")
indel_exclusion_grl <- rowRanges(indel_exclusion_vcf)
indel_exclusion_grl <- indel_exclusion_grl[seqnames(indel_exclusion_grl) %in% c(1:22, "X", "Y"),]
seqlevels(indel_exclusion_grl) <- seqlevels(indel_exclusion_grl)[1:24] # The decoy chromosomes give errors later on
seqlevels(indel_exclusion_grl) <- paste("chr", seqlevels(indel_exclusion_grl), sep = "")
indel_exclusion_grl <- get_mut_type(indel_exclusion_grl, type = "indel") # there appears to be one SNVs according to MutationalPatterns


# Filter the indels for recurrency and homopolymers
indel_grl_filter <- filter_indels(indel_grl = Indels_raw_autosomal, indel_recurrency_grl = indel_exclusion_grl, remove_homopolymer = TRUE)
indel_grl_filter <- lapply(indel_grl_filter, function(x) x[which(x$FILTER == "PASS"),])

indel_grl_flag<- filter_indels(indel_grl = Indels_raw_autosomal, indel_recurrency_grl = indel_exclusion_grl, remove_homopolymer = TRUE, mode = "flag")
indel_grl_filter <- lapply(indel_grl_filter, function(x) x[which(x$FILTER == "PASS"),])

### Count the number of SBSs in the FA samples
Indel_Counts <- data.frame(Sample = names(indel_grl_filter), Indel_RAW = lengths(indel_grl_filter))
Indel_Counts <- merge(Indel_Counts, Metadata_FA[,c("Sample", "Individual", "Age", "CallableLoci_autosomal")])
Indel_Counts$Indel_NORM <- Indel_Counts$Indel_RAW / as.numeric(Indel_Counts$CallableLoci_autosomal)

Indel_Counts$Phenotype <- "Fanconi Anemia"
Indel_Counts$Type <- ifelse(grepl("PTA", Indel_Counts$Sample) == TRUE, "PTA", "Clone")
Indel_Counts$Type <- ifelse(grepl("BULK", Indel_Counts$Sample) == TRUE, "Bulk (AML)", Indel_Counts$Type)

### Read the somatic VCFs of the HSPCs of healthy donors (which will be plotted on the ageline)
Indels_Ageline <- list()
for(Sample in Metadata_Ageline$Sample){
  print(Sample)
  vcf_file_sample <- list.files(paste(Input_dir, Metadata_Ageline$Individual[Metadata_Ageline$Sample == Sample], sep = ""), pattern =  paste(Sample,".SMuRF.filtered.vcf$", sep = ""), recursive = T, full.names = T)
  print(vcf_file_sample)
  if(length(vcf_file_sample) > 0){
    vcf_sample <- readVcf(vcf_file_sample, genome = "hg38")
    # Only take the autosomes into account
    vcf_sample <- vcf_sample[as.vector(seqnames(vcf_sample)) %in% c(1:22),]
    #SBSs_sample <- vcf_sample[isSNV(vcf_sample),]
    Indels_sample <- vcf_sample[isIndel(vcf_sample),]
    
    #SBSs_Ageline[[Sample]] <- VCF_to_GR(SBSs_sample)
    #SBSs_Ageline[[Sample]]$VAF <- as.numeric(geno(SBSs_sample)$VAF)
    
    Indels_Ageline[[Sample]] <- VCF_to_GR(Indels_sample)
    Indels_Ageline[[Sample]]$VAF <- as.numeric(geno(Indels_sample)$VAF)
  } else {
    print("!!! VCF NOT FOUND")
  }
}

Indels_Ageline <- GRangesList(Indels_Ageline)
Indels_Ageline <- get_mut_type(Indels_Ageline, type = "indel")
Indels_Ageline <- get_indel_context(Indels_Ageline, ref_genome)

Indels_all <- c(lapply(indel_grl_filter, function(x) { x[,1:8] }), as.list(Indels_Ageline))


# Count the number of SBS and indels in the healthy ageline samples
Ageline_Indels_Counts <- data.frame(Sample = names(Indels_Ageline), Indel_RAW = lengths(Indels_Ageline))
Ageline_Indels_Counts <- merge(Ageline_Indels_Counts, Metadata_Ageline, by = "Sample")

# Extrapolate the indel counts for the fraction of the genome that is not-CALLABLE
Ageline_Indels_Counts$Indel_NORM <- Ageline_Indels_Counts$Indel_RAW / as.numeric(Ageline_Indels_Counts$CallableLoci_autosomal)

Ageline_data <- Ageline_Indels_Counts[,c("Sample", "Individual", "Age", "Type", "Phenotype", "Indel_RAW", "Indel_NORM")]
Ageline_FA_Indel_Counts <- rbind(Ageline_data, Indel_Counts[,c("Sample", "Individual", "Age", "Type", "Phenotype","Indel_RAW", "Indel_NORM")])
Ageline_FA_Indel_Counts$Individual[Ageline_FA_Indel_Counts$Phenotype == "Healthy"] <- "Healthy"

# Single base substitutions ageline
lme_age = lme(Indel_NORM ~  Age, random = ~ - 1 + Age | Individual, data = Ageline_Indels_Counts)
lme_age2 = lmer(Indel_NORM ~ Age +( -1 + Age | Individual), data = Ageline_Indels_Counts)
pi = ggpredict(lme_age2,"Age",type="re")
healthy_intcpt <- lme_age$coefficients$fixed[1]
healthy_slp <- lme_age$coefficients$fixed[2]
ci_interval <- ggpredict(lme_age, ci.lvl = 0.95)$Age
head(ci_interval)
age_pval = as.data.frame(summary(lme_age)$tTable["Age","p-value"])
colnames(age_pval) <- "pval"
age_confint = intervals(lme_age)$fixed["Age",]
age_confint <- as.data.frame(t(age_confint))
age_confint$Tissue = factor("Blood")

# create data.frame with linear fits of fixed effect
indel_expected <- expand.grid(Tissue = "Blood", Age = c(min(Ageline_Indels_Counts$Age), max(Ageline_Indels_Counts$Age)))
indel_expected$fit <- predict(lme_age, level=0, newdata=indel_expected)

## These are the expected number of mutations for each patient sample (used later):
indel_expected_FA <- Metadata_FA[,c("Sample", "Age")]
indel_expected_FA$fit <- predict(lme_age, level=0, newdata=indel_expected_FA)

indel_expected_HealthyHSPCs <- Metadata_Ageline[,c("Sample", "Age")]
indel_expected_HealthyHSPCs$fit <- predict(lme_age, level=0, newdata=indel_expected_HealthyHSPCs)

Ageline_FA_Indel_Counts$fit <- predict(lme_age, level=0, newdata=Ageline_FA_Indel_Counts[,c("Sample", "Age")])


shapes_ageline <- c(23, 22, 21)
names(shapes_ageline) <- unique(Ageline_FA_Indel_Counts$Type)

sizes_ageline <- c(1.2, 1.8)
names(sizes_ageline) <- unique(Ageline_FA_Indel_Counts$Phenotype)

Ageline_FA_Indel_Counts <- merge(Ageline_FA_Indel_Counts, Metadata_FA[,c("Sample", "Mutation")], by = "Sample", all.x = T)
Ageline_FA_Indel_Counts$Mutation[is.na(Ageline_FA_Indel_Counts$Mutation)] <- ""
Ageline_FA_Indel_Counts$Individual_label <- factor(paste(Ageline_FA_Indel_Counts$Individual, " (", Ageline_FA_Indel_Counts$Mutation, ")", sep = ""))

colors_ageline <- c("darkgray", brewer.pal(n=length(unique(Ageline_FA_Indel_Counts$Individual))-1,"Dark2")) 
names(colors_ageline) <- unique(Ageline_FA_Indel_Counts$Individual)

colors_ageline2 <- c("darkgray", brewer.pal(n=length(unique(Ageline_FA_Indel_Counts$Individual))-1,"Dark2")) 
names(colors_ageline2) <- unique(Ageline_FA_Indel_Counts$Individual_label)


Indel_age_plot = ggplot() +
  geom_ribbon(data=pi, aes(x=x, ymin=conf.low, ymax=conf.high),fill = "gray90") +
  geom_ribbon(data=ci_interval, aes(x=x, ymin=conf.low, ymax=conf.high),fill = "gray60") +
  geom_line(data=indel_expected, aes(y = fit, x = Age), linewidth=0.3, col="black") +
  geom_abline(intercept = healthy_intcpt, slope = healthy_slp, color="black", linetype = "dotted", linewidth=0.3) +
  geom_point(data = Ageline_FA_Indel_Counts, aes( y = Indel_NORM, x = Age, fill = Individual_label, shape = Type, size = Phenotype), color = "black", stroke = 0.2) +
  #geom_text(data = age_pval, aes(x = 35, y = 1400,label=paste("P =",format(pval,scientific = F,digits = 2))), size=3.5, colour="black") +
  scale_y_continuous(expand = c(0,0),breaks = c(0,10,20, 30,40, 50)) +  
  scale_x_continuous(expand = c(0,0), breaks = c(0,5, 10, 15, 20))  +
  scale_fill_manual(values= colors_ageline2) +
  scale_size_manual(values = sizes_ageline, guide = "none") + 
  scale_shape_manual(values = shapes_ageline) +
  coord_cartesian(xlim = c(-2,20), ylim = c(-5,50)) +
  geom_segment(aes(x = 0, xend = 20, y = -Inf, yend = -Inf), col = "black", linewidth = 0.5) +
  geom_segment(aes(y = 0, yend = 60, x = -Inf, xend = -Inf), col = "black") +
  theme_classic(base_size = 6) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(strip.background = element_blank()) +
  theme(legend.position="right",axis.line = element_blank()) +
  labs(x = "Age (years)", y = "Indels per genome")  +
  guides(fill = guide_legend("Individual", override.aes = list(shape = 21, size = 2), keyheight = 0.5, keywidth= 0.5,  byrow = TRUE),
         shape = guide_legend("Method", override.aes = list(size = 2), keyheight = 0.5, keywidth= 0.5, byrow = TRUE)) +
  theme(legend.position="right",axis.line = element_blank(), 
        legend.spacing.y = unit(0, 'cm'), legend.margin = margin(t = 5, b = 5, l= 0, r =0, unit = "pt"),
        legend.key = element_rect(color = NA, fill = NA),
        legend.text = element_text( margin = margin(t = 0, b = 0, unit = "pt")))

Indel_age_plot

#dir.create(paste(Output_dir, "Fig3", sep = "/"))
ggsave(paste(Output_dir, "/4A_Fanconi_Ageline_Indel_zoom.pdf", sep = ""), width = 7, height = 4, units = "cm",dpi = 300)



### Observed vs expected number of single base substitutions
Indel_obs_exp <- merge(indel_expected_FA, Indel_Counts[,c("Sample", "Individual", "Indel_NORM", "Phenotype", "Type")], by = "Sample")
Indel_obs_exp$obs_exp <- Indel_obs_exp$Indel_NORM / Indel_obs_exp$fit

Indel_obs_exp_healthy <- merge(indel_expected_HealthyHSPCs, Ageline_Indels_Counts[,c("Sample", "Individual", "Indel_NORM", "Phenotype", "Type")], by = "Sample")
Indel_obs_exp_healthy$obs_exp <- Indel_obs_exp_healthy$Indel_NORM / Indel_obs_exp_healthy$fit

# only select the healthy individuals with similar ages (6-20) as the FA patients
Indel_obs_exp <- rbind(Indel_obs_exp, Indel_obs_exp_healthy[which(Indel_obs_exp_healthy$Age > 6 & Indel_obs_exp_healthy$Age < 20),])
Indel_obs_exp$Label <- Indel_obs_exp$Individual
Indel_obs_exp$Label[Indel_obs_exp$Phenotype == "Healthy"] <- "Healthy"

Indel_obs_exp$Individual <- factor(Indel_obs_exp$Individual, levels = unique(Indel_obs_exp$Individual[order(Indel_obs_exp$Age)]))
Indel_obs_exp$Label <- factor(Indel_obs_exp$Label, levels = c("Healthy", unique(Metadata_FA$Individual[order(Metadata_FA$Age)])))

wilcox <- compare_means(obs_exp ~ Label,  data = Indel_obs_exp, p.adjust.method = "bonferroni", ref.group = "Healthy")
comparisons_with_healthy <- wilcox[wilcox$group1 == "Healthy",]
# Only show the p-values for the significant differences between healthy and patients
my_comparisons <- list(c("Healthy", "PMCFANC02"), c("Healthy", "IBFM35") )

stat.test <- wilcox_test(obs_exp ~ Label, data = Indel_obs_exp, ref.group = "Healthy") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "Label", dodge = 1, step.increase = .15)
stat.test$Label <- stat.test$group2
# Only plot the significant p.adjusted values
stat.test_sig <- stat.test[stat.test$p.adj < 0.05,]
stat.test_sig$y.position[1] <- stat.test_sig$y.position[1] - 0.8
stat.test_sig$y.position[2] <- stat.test_sig$y.position[2] - 0.6

ggplot(Indel_obs_exp, aes(x = Label, y = obs_exp, fill = Label)) +
  geom_boxplot(linewidth = 0.2, alpha = 0.5, outlier.shape = NA, color = "gray30") + 
  geom_point(size = 0.5) + 
  geom_hline(yintercept = 1, linetype = 2, color = "darkgray", linewidth = 0.3) + 
  # stat_compare_means(aes(label=..p.adj..), method = "wilcox.test", 
  #                    comparisons = my_comparisons, size = 1.6, 
  #                    bracket.size = 0.2, p.adjust.methods = "BH", ref.group = "Healthy") +
  stat_pvalue_manual(
    stat.test_sig, label = paste("p=", "{p.adj}", sep = ""), size = 1.6, bracket.size = 0.2, vjust = -0.3, tip.length = 0.03) +
  theme_classic(base_size = 6) +
  scale_fill_manual(values = colors_ageline) + 
  coord_cartesian(ylim = c(0,3.6)) +
  labs(y = "Indel ratio\nobserved/expected", x = "Individual") +
  theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position = "none") 

ggsave(paste(Output_dir, "4B_Fanconi_Ageline_Indel_ObsExp.pdf", sep = ""), width = 4.5, height = 4.5, units = "cm",dpi = 300)

# In
indel_counts_FA <- count_indel_contexts(indel_grl_filter)
indel_counts_Ageline <- count_indel_contexts(Indels_Ageline)

indel_counts <- rbind(melt(indel_counts_FA), melt(indel_counts_Ageline))
colnames(indel_counts) <- c("IndelType", "Sample", "Counts_raw")
indel_counts <- merge(indel_counts, Metadata_Fig4, by = "Sample")

means_per_individual <- indel_counts %>% dplyr::group_by(Individual, IndelType) %>% 
  dplyr::summarise(sum = sum(Counts_raw)) %>% 
  mutate(Freq = sum/sum(sum))

means_per_individual2 <- means_per_individual %>% tidyr::separate(IndelType, c("muttype", "muttype_sub"), 
                                                                 sep = "_(?=[0-9])") %>% dplyr::mutate(muttype = factor(muttype, 
                                                                                                                        levels = unique(muttype)))

means_per_individual3 <- means_per_individual2  %>% dplyr::select(-muttype_sub) %>%  dplyr::group_by(Individual,muttype) %>% dplyr::summarise_all(list(~sum(.)))
means_per_individual3 <- merge(means_per_individual3, Metadata_Fig4[!duplicated(Metadata_Fig4$Individual), c("Individual", "Age", "Phenotype", "Mutation")], by = "Individual")
means_per_individual3 <- means_per_individual3[means_per_individual3$Age > 5 & means_per_individual3$Age < 20,]

indel_colors <- c("#FDBE6F", "#FF8001", "#B0DD8B", "#36A12E", "#FDCAB5", 
                  "#FC8A6A", "#F14432", "#BC141A", "#D0E1F2", "#94C4DF", 
                  "#4A98C9", "#1764AB", "#E2E2EF", "#B6B6D8", "#8683BD", 
                  "#61409B")

#counts_main$label <- ifelse(counts_main$count > max_value_to_show, counts_main$count, "")
levels(means_per_individual3$muttype) <- gsub(levels(means_per_individual3$muttype), pattern = "deletion", replacement = "del")
levels(means_per_individual3$muttype) <- gsub(levels(means_per_individual3$muttype), pattern = "insertion", replacement = "ins")
levels(means_per_individual3$muttype) <- gsub(levels(means_per_individual3$muttype), pattern = "microhomology", replacement = "MH")
means_per_individual3$Individual <- factor(means_per_individual3$Individual, levels = unique(means_per_individual3$Individual[order(means_per_individual3$Age, decreasing = F)]))

means_per_individual3$Phenotype <- factor(means_per_individual3$Phenotype, levels = c("FanconiAnemia", "AML", "Healthy"))
sums_individual <- means_per_individual3 %>%  dplyr::group_by(Individual, Phenotype) %>% dplyr::summarise(sum = sum(sum))


ggplot(means_per_individual3, aes(x = Individual, y = Freq, fill = muttype)) + 
  geom_bar(stat = "identity", col = "black", width = 0.75, linewidth = 0.2) + 
  geom_text(data = sums_individual,
            aes(y = 1.05, label = sum, x = Individual), size = 1.6, inherit.aes = F) +
  facet_grid(.~Phenotype, scales = "free", space = "free") +
  geom_segment(x = "PMCFANC08",xend = "PMCFANC02", y = 1.12,yend = 1.12, linewidth = 0.2,inherit.aes=FALSE) +
  geom_segment(x = "HSCT1",xend = "HSCT2", y = 1.12,yend = 1.12, linewidth = 0.2,inherit.aes=FALSE) +
  geom_segment(x = "IBFM35",xend = "IBFM35", y = 1.12,yend = 1.12, linewidth = 0.2,inherit.aes=FALSE) +
  #geom_text( size = 2, position = position_fill(vjust = 0.5)) +
  labs(x = "Individual", y = "Indels (fraction)", fill = "Indel type") + 
  scale_fill_manual(values = indel_colors) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.13), breaks = c(0,0.25,0.5,0.75,1)) +   
  theme_classic(base_size = 6) + 
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.2, 'cm'), 
        legend.margin = margin(t = 5, b = 5, l= 0, r =0, unit = "pt"),
        strip.background = element_blank(),strip.placement = "inside") +
  guides(fill=guide_legend(ncol=2))

ggsave(paste(Output_dir, "4C_Fanconi_Ageline_IndelSpectra.pdf", sep = ""), width = 9, height = 4, units = "cm",dpi = 300)


### 
means_per_group <- means_per_individual3
means_per_group$Group <- "FANCc-/-\nFANCA-/-"
means_per_group$Group[means_per_group$Phenotype == "Healthy"] <- "Healthy"
means_per_group$Group[means_per_group$Mutation == "BRCA2"] <- "BRCA2-/-"
means_per_group <- means_per_group[means_per_group$Phenotype != "AML",]

indel_counts %>% dplyr::group_by(Individual, IndelType) %>% 
  dplyr::summarise(sum = sum(Counts_raw)) %>% 
  mutate(Freq = sum/sum(sum))


means_per_individual2 <-merge(means_per_individual, Metadata_Fig4[!duplicated(Metadata_Fig4$Individual), c("Individual", "Age", "Phenotype", "Mutation")], by = "Individual")
means_per_individual2 <- means_per_individual2[which(means_per_individual2$Age > 5 & means_per_individual2$Age < 20),]

means_per_individual2$Group <- "FANCC-/-\nFANCA-/-"
means_per_individual2$Group[means_per_individual2$Phenotype == "Healthy"] <- "Healthy"
means_per_individual2$Group[means_per_individual2$Mutation == "BRCA2"] <- "BRCA2-/-"
means_per_group <- means_per_individual2[means_per_individual2$Phenotype != "AML",] # Exclude IBFM35 as we only have one individual
means_per_group <- means_per_group %>% dplyr::group_by(Group, IndelType) %>% 
  dplyr::summarise(meanFreq = mean(Freq)) 

indel_counts_groups <- means_per_group %>% tidyr::pivot_wider(names_from = "Group", values_from = "meanFreq") %>% tibble::column_to_rownames("IndelType")
dev.off()
pdf(file =paste(Output_dir, "4D_Fanconi_IndelContext_Grouped.pdf", sep = ""), width = 5, height = 2.5, onefile = F)
plot_indel_contexts2(indel_counts_groups[,c("Healthy", "FANCC-/-\nFANCA-/-", "BRCA2-/-")], condensed = TRUE, show_indel_counts_strips = FALSE)
dev.off()

### Figure S9C

indel_counts <- data.frame(Sample = names(indel_grl_filter), 
                           PASS = lengths(indel_grl_filter), 
                           RECURRENT = lengths(lapply(indel_grl_flag, function(x) x[which(x$FILTER =="RECURRENT"),]) ),
                           HOMOPOLYMER = lengths(lapply(indel_grl_flag, function(x) x[which(x$FILTER =="HOMOPOLYMER"),]) ),
                           RECURRENT_HOMOPOLYMER = lengths(lapply(indel_grl_flag, function(x) x[which(x$FILTER =="RECURRENT;HOMOPOLYMER"),]) ))

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
  scale_y_continuous(expand = c(0,0), limit= c(0,900), breaks = seq(0,900, 300)) +
  theme_classic(base_size = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.25, 'cm')) +
  labs(x = "Sample", y = "Indels (VAF>0.25)", fill = "Filter")

ggsave(paste(Output_dir_FigureS9, "FigureS9C_Indel_Filter.pdf", sep = ""), width = 14, height = 5, units = "cm")



## Deletion sizes (figure 5)


SVs <- read.delim("~/hpc/pmc_vanboxtel/projects/PTATO/3_Output/Somatic_Data/Table_S3.txt")

SV_list <- list()

for(Sample in unique(SVs$Sample)){
  print(Sample)
  SVs_sample <- GRanges(seqnames = paste("chr", SVs$Chr1[SVs$Sample == Sample], sep = ""), IRanges(start = SVs$Start1[SVs$Sample == Sample],
                                                                                                   end = SVs$End1[SVs$Sample == Sample]))
  # Some deletions detected by the SV pipeline may already be called by GATK. Remove these from the SV list
  Indels_sample <- Indels_all[[Sample]]
  if(!is.null(Indels_sample)){
    Overlap_Indels_SVs <- findOverlaps(SVs_sample, Indels_sample)
    if(length(Overlap_Indels_SVs) > 0){
      SVs_sample <- SVs_sample[-queryHits(Overlap_Indels_SVs)]
    }
    SVs_sample <- SVs_sample[seqnames(SVs_sample) %in% paste("chr", 1:22, sep = "")]
    SVs_sample$size <- width(SVs_sample)
    SV_list[[Sample]] <- SVs_sample
  }
}


Indels_all[["PMCFANC02-BMHSPC-1C22"]][seqnames(Indels_all[["PMCFANC02-BMHSPC-1C22"]]) == "chr22",]
Indels_all[["PMCFANC02-BMHSPC-1L5"]][seqnames(Indels_all[["PMCFANC02-BMHSPC-1L5"]]) == "chrX",]

indel_sizes <- data.frame()

# collect the size of each indel per sample
for(individual in unique(Metadata_Fig4$Individual)){
  print(individual)
  for(sample in Metadata_Fig4$Sample[Metadata_Fig4$Individual == individual]){
    print(sample)
    indels_sample <- Indels_all[[sample]]
    if(!is.null(indels_sample)){
      indel_size_sample <- data.frame(Individual = individual, Sample = sample, 
                                      Type = ifelse(lengths(indels_sample$REF) > lengths(unlist(indels_sample$ALT)), "DEL", "INS"),
                                      size = abs(lengths(indels_sample$REF) - lengths(unlist(indels_sample$ALT))))
      
      if(!is.null(SV_list[[sample]])){
        deletion_size_sample <- data.frame(Individual = individual, Sample = sample, 
                                        Type = "DEL",
                                        size = SV_list[[sample]]$size)
        indel_size_sample <- rbind(indel_size_sample, deletion_size_sample)
      }
      
      if(nrow(indel_sizes) > 0 ){
        indel_sizes <- rbind(indel_sizes, indel_size_sample)
      } else {
        indel_sizes <- indel_size_sample
      }
    }

  }
}

indel_sizes$Category <- ifelse(indel_sizes$size > 1, "2-10bp", "1bp")
indel_sizes$Category <- ifelse(indel_sizes$size > 10, "11-50bp", indel_sizes$Category)
indel_sizes$Category <- ifelse(indel_sizes$size > 50, ">50bp", indel_sizes$Category)
indel_sizes$Category <- factor(indel_sizes$Category, levels = c(">50bp", "11-50bp","2-10bp", "1bp"))

# PEr sample
indel_sizes_PerSample <- indel_sizes[indel_sizes$Type == "DEL",] %>% group_by(Sample, Category, .drop = F) %>% 
  summarise(Freq = n(), counts = n()) %>%
  mutate(Freq = Freq/sum(Freq))
indel_sizes_PerSample <- merge(indel_sizes_PerSample, Metadata_Fig4[,c("Sample", "Individual", "Age", "Phenotype")], by = "Sample", all.x = T)

# only show the age-matched samples
indel_sizes_PerSample <- indel_sizes_PerSample[indel_sizes_PerSample$Age > 6 & indel_sizes_PerSample$Age < 20,]

# mean per individual
indel_sizes_perIndividual <- indel_sizes_PerSample %>% 
  group_by(Individual, Category) %>% 
  summarise(mean = mean(Freq), sum = sum(counts), samples = length(unique(Sample)))
indel_sizes_perIndividual <- merge(indel_sizes_perIndividual, 
                                   Metadata_Fig4[!duplicated(Metadata_Fig4$Individual),c("Individual", "Age", "Phenotype")], by = "Individual", all.x = T)

# Per phenotype

indel_sizes_perGroup <- as.data.frame(indel_sizes_perIndividual %>% 
                                        group_by(Phenotype, Category)  %>% 
                                        summarise(mean = mean(mean), sum = sum(sum), samples = sum(samples)))
indel_sizes_perGroup$Species <- "Human"
# This data is obtained from Garaycoechea et al, Nature (2018). Figure 5f (supplemental table)
indel_sizes_Mice <- data.frame(Phenotype = c(rep("Healthy", 4),
                                                 rep("Fancd2-/-",4),
                                             rep("Aldh2-/-;\nFancd2-/-", 4)), 
                               Category = rep(c(">50bp", "11-50bp", "2-10bp", "1bp"),3),
                               sum = c(0, 2, 1, 6,
                                       0, 10, 6, 19,
                                       17,30,18,34),
                               mean = c(0, 0.222, 0.111, 0.667,
                                 0, 0.286, 0.171, 0.543,
                                        0.172, 0.303, 0.182, 0.343),
                               samples = c(rep(3,4), rep(4,4), rep(5,4)),
                               Species = "Mouse")

indel_sizes_FA <- rbind(indel_sizes_perGroup, indel_sizes_Mice)
indel_sizes_FA$Phenotype <- factor(indel_sizes_FA$Phenotype, levels = c("Healthy", "FanconiAnemia", "AML", "Fancd2-/-", "Aldh2-/-;\nFancd2-/-"))
indel_sizes_FA <- indel_sizes_FA[indel_sizes_FA$Phenotype != "AML",]
indel_sizes_FA$Label <- paste(indel_sizes_FA$Phenotype, "\n" , "(n=",indel_sizes_FA$samples, ")", sep = "")
indel_sizes_FA$Label <- factor(indel_sizes_FA$Label, levels = unique(indel_sizes_FA$Label[order(indel_sizes_FA$Phenotype)]))

total_indel_counts <- indel_sizes_FA %>% 
  group_by(Label, Species)  %>% 
  summarise( sum = sum(sum))

ggplot(indel_sizes_FA, aes(x = Label, y = mean)) + 
  geom_bar(stat = "identity", aes(fill = Category), col = "black", size = 0.2, width = 0.8) +
  facet_grid(.~Species, scales = "free", space = "free") +
  geom_text(data = total_indel_counts, aes(label = sum, y = 1.05, x = Label), size = 1.6) +
  geom_segment(x =levels(indel_sizes_FA$Label)[1],xend = levels(indel_sizes_FA$Label)[3], y = 1.1,yend = 1.1, size = 0.2,inherit.aes=FALSE) +
  geom_segment(x = levels(indel_sizes_FA$Label)[2],xend = levels(indel_sizes_FA$Label)[5], y = 1.1,yend = 1.1, size = 0.2,inherit.aes=FALSE) + 
  theme_classic(base_size  = 6) + 
  scale_fill_brewer(palette = "Reds", direction = -1) +
  scale_y_continuous(limits = c(0, 1.1), expand = c(0,0)) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.key.size = unit(0.2, 'cm'),
        legend.position = "bottom", strip.background = element_blank()) +
  labs(y = "Relative contribution deletions", fill = "Deletion size", x = "Type of HSPC")
ggsave(paste(Output_dir, "5h_Fanconi_IndelSizes.pdf", sep = ""), width = 4.5, height = 5, units = "cm",dpi = 300)
