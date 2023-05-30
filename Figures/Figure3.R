### Figures 3 and S9A-S9B

library(MutationalPatterns)
library(ggplot2)
library(VariantAnnotation)
library(RColorBrewer)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
library(nlme)
library(ggeffects)
library(lme4)
library(ggpubr)
library(reshape2)
library(ggh4x)
library(rstatix)

Input_dir <- "/path/to/MendeleyData/directory/"
Output_dir <-  "/path/to/output/directory/"
Output_dir_fig3 <- paste(Output_dir, "Figure3/", sep = "")
if(dir.exists(Output_dir_fig3) == F){ dir.create(Output_dir_fig3)}
Scripts_dir <- "/path/to/GithubScripts/directory/"

source(paste(Scripts_dir, "GeneralFunctions.R", sep = ""))

Metadata <- read.delim(paste(Input_dir, "Table_S1.txt", sep = ""))
Metadata_FA <- Metadata[Metadata$Phenotype == "FanconiAnemia" | Metadata$Individual == "IBFM35",]
Metadata_FA <- Metadata_FA[Metadata_FA$Germline == FALSE,]
Metadata_Ageline <- Metadata[Metadata$Phenotype == "Healthy" & Metadata$Training == "No" & Metadata$Type == "Clone" & Metadata$Tissue == "BoneMarrow",]
Metadata_Fig3 <- rbind(Metadata_FA,Metadata_Ageline)

### Read the PTATO vcfs
# This version of PTATO contains an error calculating the final PTAprobsCutoff. Recalculate
ptatotables_files <- list.files(Input_dir, pattern = ".ptatotable.txt", recursive = T, full.names = T)
new_cutoffs <- data.frame()
for(ptatotable_file in ptatotables_files){
  print(ptatotable_file)
  ptatotable <- read.table(ptatotable_file, header = T)
  ptatotable$prec_recall <- abs(ptatotable$precision-ptatotable$recall)
  
  ptatotable <- ptatotable[ptatotable$precision > 0 & ptatotable$recall > 0,]
  
  lowest <- ptatotable[which(ptatotable$prec_recall == min(ptatotable$prec_recall, na.rm = T)),]
  if(nrow(lowest) > 1){
    lowest <- lowest[1,]
  }
  
  new_cutoff_sample <- data.frame(file = ptatotable_file, walker_cutoff = lowest$PTAprob_cutoff)
  new_cutoffs <- rbind(new_cutoffs, new_cutoff_sample)
}


SBSs_raw <- list()
Indels_raw <- list()
MinimalVAF <- 0.2 # all variants will have a VAF higher than MinimalVAF
for(Sample in Metadata_FA$Sample){
  print(Sample)
  
  # For PTA samples we use the PTATO-filtered VCFs, for bulk samples we use the SMuRF-filtered VCFs
  if(Metadata_FA$Type[Metadata_FA$Sample == Sample] == "PTA"){
    # PTATO SBSs
    vcf_unfiltered <- list.files(paste(Input_dir, Metadata_FA$Individual[Metadata_FA$Sample == Sample], "/PTATO/", sep = ""),
                                 pattern = paste(Sample, ".snvs.ptato.callable.vcf.gz$", sep = ""), full.names = T, recursive = T)
    print(vcf_unfiltered)
    vcf_filtered <- list.files(paste(Input_dir, Metadata_FA$Individual[Metadata_FA$Sample == Sample], "/PTATO/", sep = ""),
                                 pattern = paste(Sample, ".snvs.ptato.filtered.vcf.gz$", sep = ""), full.names = T, recursive = T)
    print(vcf_filtered)
    SBSs_raw[[Sample]] <- readPTATOvcf(vcf_unfiltered = vcf_unfiltered, 
                                      vcf_filtered = vcf_filtered, 
                                      VAF_threshold = MinimalVAF)
    
    SBSs_raw[[Sample]]$PTAprobCutoff_Walker <- new_cutoffs$walker_cutoff[grep(pattern = Sample, x= new_cutoffs$file)]
    SBSs_raw[[Sample]]$PTAprobCutoff <- (unique(as.numeric(SBSs_raw[[Sample]]$PTAprobCutoff_Walker)) + 
                                                  unique(as.numeric(SBSs_raw[[Sample]]$PTAprobCutoff_Cossim))) / 2
    SBSs_raw[[Sample]]$FILTER <- ifelse(SBSs_raw[[Sample]]$PTAprob > SBSs_raw[[Sample]]$PTAprobCutoff, "FAIL", "PASS")
    SBSs_raw[[Sample]]$FILTER[ SBSs_raw[[Sample]]$VAF < 0.2] <- "FAIL_VAF"
    
    # PTATO Indels
    indels_vcf_file <- list.files(paste(Input_dir, Metadata_FA$Individual[Metadata_FA$Sample == Sample], "/PTATO/", sep = ""),
                             pattern = paste(Sample, ".indels.ptato.filtered.vcf.gz$", sep = ""), full.names = T, recursive = T)
    print(indels_vcf_file)
    Indels_Sample <- readVcf(indels_vcf_file, genome = "hg38")
    Indels_Sample <- Indels_Sample[as.vector(seqnames(Indels_Sample)) %in% c(1:22, "X", "Y"),]
    
    Indels_raw[[Sample]] <- VCF_to_GR(Indels_sample)
    Indels_raw[[Sample]]$VAF <- geno(Indels_sample)$VAF
    
  } else if(Metadata_FA$Type[Metadata_FA$Sample == Sample]%in% c("Clone","Bulk")){
    # SMuRF VCFs contain both somatic SBSs and indels
    vcf_file <- list.files(paste(Input_dir, Metadata_FA$Individual[Metadata_FA$Sample == Sample], "/PTATO/", sep = ""),
                      pattern = paste(Sample, ".SMuRF.filtered.callable.vcf.gz$", sep = ""), full.names = T, recursive = T)
    print(vcf_file)
    vcf <- readVcf(vcf_file, genome = "hg38")
    vcf <- vcf[as.vector(seqnames(vcf)) %in% c(1:22, "X", "Y"),]
    SBSs_sample <- vcf[isSNV(vcf),]
    Indels_sample <- vcf[isIndel(vcf),]
    
    SBSs_raw[[Sample]] <- VCF_to_GR(SBSs_sample)
    SBSs_raw[[Sample]]$VAF <- geno(SBSs_sample)$VAF
    Indels_raw[[Sample]] <- VCF_to_GR(Indels_sample)
    Indels_raw[[Sample]]$VAF <- geno(Indels_sample)$VAF
  }
}

# Only autosomal
SBSs_raw_autosomal <-  lapply(SBSs_raw, function(x) x[which(as.vector(seqnames(x)) %in% paste("chr", c(1:22), sep = "")),])

# Remove all variants that are UNCALLABLE and/or with a PTAprobs < PTAprobsCutoff
SBSs_PASS <- lapply(SBSs_raw_autosomal, function(x) x[which(x$FILTER =="PASS"),])
#sSNVs_FAIL <- lapply(SBSs_raw_autosomal, function(x) x[which(x$FILTER =="FAIL"),])
SBSs_FAIL <- lapply(SBSs_raw_autosomal, function(x) x[which(x$FILTER =="FAIL"),])

cutoffs <- lapply(SBSs_raw_autosomal, function(x) unique(x$PTAprobCutoff))


### Count the number of SBSs in the FA samples
Variant_Counts <- Metadata_FA
Variant_Counts <- data.frame(Sample = names(SBSs_PASS), SBS_RAW = lengths(SBSs_PASS), SBS_FAIL = lengths(SBSs_FAIL))
Variant_Counts <- merge(Variant_Counts, Metadata_FA[,c("Sample", "Label", "Label2", "Individual", "Age", "CallableLoci_autosomal")])
Variant_Counts$SBS_NORM <- Variant_Counts$SBS_RAW / as.numeric(Variant_Counts$CallableLoci_autosomal)
Variant_Counts$Phenotype <- "FanconiAnemia"
Variant_Counts$Type <- ifelse(grepl("PTA", Variant_Counts$Sample) == TRUE, "PTA", "Clone")
Variant_Counts$Type <- ifelse(grepl("BULK", Variant_Counts$Sample) == TRUE, "Bulk (AML)", Variant_Counts$Type)


### Read the somatic VCFs of the HSPCs of healthy donors (which will be plotted on the ageline)
SBSs_Ageline <- list()
Indels_Ageline <- list()
for(Sample in Metadata_Ageline$Sample){
  print(Sample)
  vcf_file_sample <- list.files(paste(Input_dir, Metadata_Ageline$Individual[Metadata_Ageline$Sample == Sample], sep = ""), pattern =  paste(Sample,".SMuRF.filtered.vcf$", sep = ""), recursive = T, full.names = T)
  print(vcf_file_sample)
  if(length(vcf_file_sample) > 0){
    vcf_sample <- readVcf(vcf_file_sample, genome = "hg38")
    # Only take the autosomes into account
    vcf_sample <- vcf_sample[as.vector(seqnames(vcf_sample)) %in% c(1:22),]
    SBSs_sample <- vcf_sample[isSNV(vcf_sample),]
    Indels_sample <- vcf_sample[isIndel(vcf_sample),]
    
    SBSs_Ageline[[Sample]] <- VCF_to_GR(SBSs_sample)
    SBSs_Ageline[[Sample]]$VAF <- as.numeric(geno(SBSs_sample)$VAF)
    
    Indels_Ageline[[Sample]] <- VCF_to_GR(Indels_sample)
    Indels_Ageline[[Sample]]$VAF <- as.numeric(geno(Indels_sample)$VAF)
  } else {
    print("!!! VCF NOT FOUND")
  }

}

SBSs_all <- c(SBSs_PASS, SBSs_Ageline)

# Count the number of SBS and indels in the healthy ageline samples
Ageline_SBS_Counts <- data.frame(Sample = names(SBSs_Ageline), SBS_RAW = lengths(SBSs_Ageline), Indel_Raw = lengths(Indels_Ageline))
Ageline_SBS_Counts <- merge(Ageline_SBS_Counts, Metadata_Ageline, by = "Sample")

# Extrapolate the SBS and indel counts for the fraction of the genome that is not-CALLABLE
Ageline_SBS_Counts$SBS_NORM <- Ageline_SBS_Counts$SBS_RAW / as.numeric(Ageline_SBS_Counts$CallableLoci_autosomal)
Ageline_SBS_Counts$INDEL_NORM <- Ageline_SBS_Counts$Indel_Raw / as.numeric(Ageline_SBS_Counts$CallableLoci_autosomal)

Ageline_data <- Ageline_SBS_Counts[,c("Sample", "Individual", "Age", "Type", "Phenotype", "SBS_RAW", "SBS_NORM")]
Ageline_SBS_Counts_FA <- rbind(Ageline_data, Variant_Counts[,c("Sample", "Individual", "Age", "Type", "Phenotype","SBS_RAW", "SBS_NORM")])
Ageline_SBS_Counts_FA$Individual[Ageline_SBS_Counts_FA$Phenotype == "Healthy"] <- "Healthy"

# Single base substitutions ageline
lme_age = lme(SBS_NORM ~  Age, random = ~ - 1 + Age | Individual, data = Ageline_SBS_Counts)
lme_age2 = lmer(SBS_NORM ~ Age +( -1 + Age | Individual), data = Ageline_SBS_Counts)
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
mut_expected <- expand.grid(Tissue = "Blood", Age = c(min(Ageline_SBS_Counts$Age), max(Ageline_SBS_Counts$Age)))
mut_expected$fit <- predict(lme_age, level=0, newdata=mut_expected)

## These are the expected number of mutations for each patient sample (used later):
mut_expected_FA <- Metadata_FA[,c("Sample", "Age")]
mut_expected_FA$fit <- predict(lme_age, level=0, newdata=mut_expected_FA)

mut_expected_HealthyHSPCs <- Metadata_Ageline[,c("Sample", "Age")]
mut_expected_HealthyHSPCs$fit <- predict(lme_age, level=0, newdata=mut_expected_HealthyHSPCs)

Variant_Counts$fit <- as.numeric(predict(lme_age, level=0, newdata=Variant_Counts[,c("Sample", "Age")]))
Ageline_SBS_Counts_FA$fit <- predict(lme_age, level=0, newdata=Ageline_SBS_Counts_FA[,c("Sample", "Age")])

Ageline_SBS_Counts_FA <- merge(Ageline_SBS_Counts_FA, Metadata_FA[,c("Sample", "Mutation")], by = "Sample", all.x = T)
Ageline_SBS_Counts_FA$Mutation[is.na(Ageline_SBS_Counts_FA$Mutation)] <- ""
Ageline_SBS_Counts_FA$Individual_label <- factor(paste(Ageline_SBS_Counts_FA$Individual, " (", Ageline_SBS_Counts_FA$Mutation, ")", sep = ""))

colors_ageline <- c("darkgray", brewer.pal(n=length(unique(Ageline_SBS_Counts_FA$Individual))-1,"Dark2")) 
names(colors_ageline) <- unique(Ageline_SBS_Counts_FA$Individual)

colors_ageline2 <- c("darkgray", brewer.pal(n=length(unique(Ageline_SBS_Counts_FA$Individual))-1,"Dark2")) 
names(colors_ageline2) <- unique(Ageline_SBS_Counts_FA$Individual_label)


shapes_ageline <- c(23, 22, 21)
names(shapes_ageline) <- unique(Ageline_SBS_Counts_FA$Type)

sizes_ageline <- c(1.2, 1.8)
names(sizes_ageline) <- unique(Ageline_SBS_Counts_FA$Phenotype)

SBS_age_plot = ggplot() +
  geom_ribbon(data=pi, aes(x=x, ymin=conf.low, ymax=conf.high),fill = "gray90") +
  geom_ribbon(data=ci_interval, aes(x=x, ymin=conf.low, ymax=conf.high),fill = "gray60") +
  geom_line(data=mut_expected, aes(y = fit, x = Age), linewidth=0.3, col="black") +
  geom_abline(intercept = healthy_intcpt, slope = healthy_slp, color="black", linetype = "dotted", linewidth=0.3) +
  geom_point(data = Ageline_SBS_Counts_FA, aes( y = SBS_NORM, x = Age, fill = Individual_label, shape = Type, size = Phenotype), color = "black", stroke = 0.2) +
  ylab("Base substitutions per genome") +
  #geom_text(data = age_pval, aes(x = 35, y = 1400,label=paste("P =",format(pval,scientific = F,digits = 2))), size=3.5, colour="black") +
  scale_y_continuous(expand = c(0,0),breaks = c(0,250, 500, 750, 1000, 1250)) +  
  scale_x_continuous(expand = c(0,0), breaks = c(0,5, 10, 15, 20))  +
  scale_fill_manual(values= colors_ageline2) +
  scale_size_manual(values = sizes_ageline, guide = "none") + 
  scale_shape_manual(values = shapes_ageline) +
  coord_cartesian(xlim = c(-2,20), ylim = c(-65,1100)) +
geom_segment(aes(x = 0, xend = 70, y = -Inf, yend = -Inf), col = "black", linewidth = 0.5) +
  geom_segment(aes(y = 0, yend = 1500, x = -Inf, xend = -Inf), col = "black") +
  theme_classic(base_size = 6) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(strip.background = element_blank()) +
  theme(legend.position="right",axis.line = element_blank()) +
  xlab("Age (years)")  +
  guides(fill = guide_legend("Individual", override.aes = list(shape = 21, size = 2), keyheight = 0.5, keywidth= 0.5,  byrow = TRUE),
         shape = guide_legend("Method", override.aes = list(size = 2), keyheight = 0.5, keywidth= 0.5, byrow = TRUE)) +
  theme(legend.position="right",axis.line = element_blank(), 
        legend.spacing.y = unit(0, 'cm'), legend.margin = margin(t = 5, b = 5, l= 0, r =0, unit = "pt"),
        legend.key = element_rect(color = NA, fill = NA),
        legend.text = element_text( margin = margin(t = 0, b = 0, unit = "pt")))

SBS_age_plot

#dir.create(paste(Output_dir, "Fig3", sep = "/"))
ggsave(paste(Output_dir_fig3, "3a_Fanconi_Ageline_SNV_zoom.pdf", sep = ""), width = 6.8, height = 4, units = "cm",dpi = 300)

### Observed vs expected number of single base substitutions
SBS_obs_exp <- merge(mut_expected_FA, Variant_Counts[,c("Sample", "Individual", "SBS_NORM", "Phenotype", "Type")], by = "Sample")
SBS_obs_exp$obs_exp <- SBS_obs_exp$SBS_NORM / SBS_obs_exp$fit

SBS_obs_exp_healthy <- merge(mut_expected_HealthyHSPCs, Ageline_SBS_Counts[,c("Sample", "Individual", "SBS_NORM", "Phenotype", "Type")], by = "Sample")
SBS_obs_exp_healthy$obs_exp <- SBS_obs_exp_healthy$SBS_NORM / SBS_obs_exp_healthy$fit

# only select the healthy individuals with similar ages (6-20) as the FA patients
SBS_obs_exp <- rbind(SBS_obs_exp, SBS_obs_exp_healthy[which(SBS_obs_exp_healthy$Age > 6 & SBS_obs_exp_healthy$Age < 20),])
SBS_obs_exp$Label <- SBS_obs_exp$Individual
SBS_obs_exp$Label[SBS_obs_exp$Phenotype == "Healthy"] <- "Healthy"

SBS_obs_exp$Individual <- factor(SBS_obs_exp$Individual, levels = unique(SBS_obs_exp$Individual[order(SBS_obs_exp$Age)]))
SBS_obs_exp$Label <- factor(SBS_obs_exp$Label, levels = c("Healthy", unique(Metadata_FA$Individual[order(Metadata_FA$Age)])))

# Showing multiple testing adjusted p-values does not work for the standard stat_compare_means function
wilcox <- compare_means(obs_exp ~ Label,  data = SBS_obs_exp)
comparisons_with_healthy <- wilcox[wilcox$group1 == "Healthy",]
# Only show the p-values for the significant differences between healthy and patients
my_comparisons <- list( c("Healthy", "PMCFANC03"), c("Healthy", "PMCFANC02"), c("Healthy", "IBFM35") )

stat.test <- wilcox_test(obs_exp ~ Label, data = SBS_obs_exp, ref.group = "Healthy") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "Label", dodge = 1, step.increase = .13)
stat.test$Label <- stat.test$group2
# Only plot the significant p.adjusted values
stat.test_sig <- stat.test[stat.test$p.adj < 0.05,]
stat.test_sig$y.position <- stat.test_sig$y.position - 1.1

ggplot(SBS_obs_exp, aes(x = Label, y = obs_exp, fill = Label)) +
  geom_boxplot(linewidth = 0.2, alpha = 0.5, outlier.shape = NA, color = "gray30") + 
  geom_point(size = 0.5) + 
  geom_hline(yintercept = 1, linetype = 2, color = "darkgray", size = 0.3) + 
  #stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, size = 1.6, bracket.size = 0.2) +
  stat_pvalue_manual(
    stat.test_sig, label = paste("p=", "{p.adj}", sep = ""), size = 1.6, bracket.size = 0.2, vjust = 0, tip.length = 0.03) +
  theme_classic(base_size = 6) +
  scale_fill_manual(values = colors_ageline) + 
  coord_cartesian(ylim = c(0,4)) +
  labs(y = "Ratio base substitutions\nobserved/expected", x = "Individual") +
  theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position = "none") 

ggsave(paste(Output_dir_fig3, "/3b_Fanconi_Ageline_sSNV_ObsExp.pdf", sep = ""), width = 4, height = 4, units = "cm",dpi = 300)

### Mutation spectra
type_occurrences <- mut_type_occurrences(lapply(SBSs_PASS, "[", , 1:5), ref_genome) # The columns are different between PTA and non-PTA samples

type_occurrences$Sample <- row.names(type_occurrences)
type_occurrences <- merge(type_occurrences, Metadata_FA[,c("Sample", "Label", "Individual", "Age")], by = "Sample")
type_occurrences$Phenotype <- "Fanconi anemia"

type_occurrences_ageline <- mut_type_occurrences(SBSs_Ageline, ref_genome)
type_occurrences_ageline$Sample <- row.names(type_occurrences_ageline)
type_occurrences_ageline <- merge(type_occurrences_ageline, Metadata_Ageline[,c("Sample", "Label", "Individual", "Age")], by = "Sample")
type_occurrences_ageline$Phenotype <- "Healthy"

# Only include the healthy donors within the same age range as the patients with FA
type_occurrences_merged <- rbind(type_occurrences, type_occurrences_ageline[type_occurrences_ageline$Age > 5 &  type_occurrences_ageline$Age < 15,])
type_occurrences_merged_cpg <- type_occurrences_merged[,c(2,3,8,9,5,6, 10, 11,12,13)]
row.names(type_occurrences_merged_cpg) <- type_occurrences_merged$Sample

type_occurrences_merged2 <- type_occurrences_merged_cpg[,c(1:6)]
by = type_occurrences_merged_cpg$Individual

# Change to long format for ggplot
tb_per_sample <-type_occurrences_merged2 %>% tibble::rownames_to_column("sample") %>% 
  dplyr::mutate(by = by) %>% tidyr::pivot_longer(c(-sample, 
                                                   -by), names_to = "variable", values_to = "nmuts") %>% 
  dplyr::group_by(sample) %>% dplyr::mutate(value = nmuts/sum(nmuts)) %>% 
  dplyr::ungroup() %>% dplyr::mutate(sub_type = stringr::str_remove(variable, 
                                                                    " .*"), variable = factor(variable, levels = unique(variable)))

tb <- tb_per_sample %>% dplyr::mutate(by = factor(by, levels = unique(by))) %>% 
  dplyr::group_by(by, variable) %>% dplyr::summarise(sub_type = sub_type[[1]], 
                                                     mean = mean(value), stdev = stats::sd(value), total_individuals = sum(value), 
                                                     total_mutations = sum(nmuts)) %>% dplyr::mutate(total_individuals = sum(total_individuals), 
                                                                                                     total_mutations = sum(total_mutations)) %>% dplyr::mutate(sem = stdev/sqrt(total_individuals), error_95 = ifelse(total_individuals > 1, qt(0.975, df = total_individuals -  1) * sem, NA)) %>% dplyr::ungroup() %>% dplyr::mutate(total_mutations = total_mutations, error_pos = mean)

colnames(tb)[1] <- "Individual"
tb2 <- merge(tb, type_occurrences_merged_cpg[!duplicated(type_occurrences_merged_cpg$Individual),c("Individual","Age","Phenotype")], by = "Individual", all.x = T)
total_muts <- tb2[!duplicated(tb2$Individual),]
# Reorder samples based on age
tb2$Individual <- factor(tb2$Individual, levels = unique(tb2$Individual[order(tb2$Age, decreasing = F)]))

ggplot(tb2, aes(x = Individual, fill = variable, y = mean)) + 
  geom_bar(stat = "identity", width = 0.8) +
  geom_text(data = total_muts, aes(label = total_mutations, y = 1.05), size = 1.5) +
  # Lines instead of rectangles as facet strips:
  geom_segment(x = "PMCFANC08",xend = "IBFM35", y = 1.12,yend = 1.12, size = 0.2,inherit.aes=FALSE) +
  geom_segment(x = "HSCT1",xend = "HSCT2", y = 1.12,yend = 1.12, size = 0.2,inherit.aes=FALSE) +
  facet_grid(.~Phenotype, scales = "free_x", space = "free_x", ) +
  scale_fill_manual(values = COLORS7) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.13), labels = scales::percent, breaks = seq(0,1,0.25)) +
  theme_classic(base_size = 6) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        legend.key.size = unit(0.25, 'cm'), strip.background = element_blank(),strip.placement = "inside") +
  labs(y = "Relative contribution", fill = "Type of base\nsubstitution")

ggsave(paste(Output_dir_fig3, "3c_Fanconi_MutationSpectrum.pdf", sep = ""), width = 7, height = 4, units = "cm",dpi = 300)

### 96-profiles
mut_mat_merged <- mut_matrix(vcf_list = lapply(SBSs_all, "[", , 1:5), ref_genome = ref_genome)

# Calculate the relative contribution of each channel to the 96-profile
mut_mat_rel <- sweep(mut_mat_merged,2,colSums(mut_mat_merged),`/`)

# Calculate the mean relative contribution per individual
mut_mat_individuals_abs <- matrix(nrow = nrow(mut_mat_merged), ncol = 0)
row.names(mut_mat_individuals_abs) <- row.names(mut_mat_individuals_abs)
mut_mat_individuals_rel <- matrix(nrow = nrow(mut_mat_rel), ncol = 0)
row.names(mut_mat_individuals_rel) <- row.names(mut_mat_individuals_rel)
for(individual in unique(Metadata_Fig3$Individual)){
  print(individual)
  mut_mat_individual_abs <- mut_mat_merged[,Metadata_Fig3$Sample[Metadata_Fig3$Individual == individual]]
  mut_mat_individual_rel <- mut_mat_rel[,Metadata_Fig3$Sample[Metadata_Fig3$Individual == individual]]
  if(!is.null(ncol(mut_mat_individual_rel))){
    mut_mat_individuals_rel <- cbind(mut_mat_individuals_rel, rowMeans(mut_mat_individual_rel))
    mut_mat_individuals_abs <- cbind(mut_mat_individuals_abs, rowSums(mut_mat_individual_abs))
    
  } else {
    mut_mat_individuals_rel <- cbind(mut_mat_individuals_rel,mut_mat_individual_rel)
    mut_mat_individuals_abs <- cbind(mut_mat_individuals_abs,mut_mat_individual_abs)
    
  }
  colnames(mut_mat_individuals_rel)[ncol(mut_mat_individuals_rel)] <- individual
  colnames(mut_mat_individuals_abs)[ncol(mut_mat_individuals_abs)] <- individual
}

# Calculate the mean relative contributions per subgroup 
Groups <- list(Healthy = c("HSCT1", "HSCT2", "HSCT3", "HSCT4"),
               `FANCC-/-\nFANCA-/-` = c("PMCFANC01", "PMCFANC03" ,  "PMCFANC06"  , "PMCFANC08"),
               `BRCA2-/-` = c("PMCFANC02"))

mut_mat_grouped <- matrix(nrow = nrow(mut_mat_rel), ncol = 0)
row.names(mut_mat_grouped) <- row.names(mut_mat_rel)

for(group in names(Groups)){
  print(group)
  mut_mat_group <- mut_mat_individuals_rel[,Groups[[group]]]
  if(!is.null(ncol(mut_mat_group))){
    mut_mat_grouped <- cbind(mut_mat_grouped, rowMeans(mut_mat_group))
  } else {
    mut_mat_grouped <- cbind(mut_mat_grouped,mut_mat_group)
  }
  colnames(mut_mat_grouped)[ncol(mut_mat_grouped)] <- group
}

dev.off()
pdf(paste(Output_dir_fig3, "3d_FA_96profiles_merged.pdf", sep = ""), width = 3, height = 2, onefile = F)
plot_96_profile_small(mut_mat_grouped, condensed = T, ymax = 0.1, add_total = FALSE) 
dev.off()

### Signatures
signatures <- get_known_signatures()
signatures <- signatures[,1:30]
pta_signature <- read.table(paste(Input_dir, "/Resources/PTA_Artefact_Signature.txt", sep = ""), header=T)
signatures <- cbind(signatures,pta_signature$PTA)
colnames(signatures)[length(colnames(signatures))] <- "PTA"
#Read in the existing signatures.
cosmic_sig_fname = paste(Input_dir, "/Resources/sigProfiler_SBS_working_signatures_incl_hspc.txt", sep = "")
cosmic_sign = read.table(cosmic_sig_fname, sep = "\t", header = T)
cosmic_sign = as.matrix(cosmic_sign[,-c(1,2)])
HSPC <- cosmic_sign[,"HSPC"]
signatures <- cbind(HSPC, signatures)
# Signatures commonly found in hematopoietic cells (eg Osorio et al, Brandsma et al, Mitchell et al)
# SBS HSPC is also known as SBSblood 
signatures_blood <- signatures[,c("HSPC", "SBS1", "SBS3", "SBS5", "SBS18", "PTA")]
# Merge signatures that are very similar
merged_signatures <- merge_signatures(signatures, cos_sim_cutoff = 0.8)

# Perform strict signature refit using 100 bootstraps (and subsequently taking the mean contribution from each bootstrap)
contri_boots <- fit_to_signatures_bootstrapped(mut_mat_merged,
                                               signatures_blood,
                                               n_boots = 100,
                                               method = "strict", max_delta = 0.005)
#plot_bootstrapped_contribution(contri_boots, mode = "relative", plot_type = "bar")

# Change to long format for plotting
contri_tb <- contri_boots %>% as.data.frame() %>% tibble::rownames_to_column("exp") %>% 
  tidyr::gather(key = "sig", value = "contri", -exp) %>% 
  dplyr::mutate(Sample = gsub("_[^_]+$", "", exp), Sample = factor(Sample, 
                                                                   levels = unique(Sample)), sig = factor(sig, levels = unique(sig)))
# Calculate the mean contribution for each signature
contri_tb2 <- contri_tb %>% dplyr::group_by(Sample, sig) %>% 
  dplyr::summarise(mean = mean(contri)) %>% 
  dplyr::ungroup()

# Calculate the number of mutations that are not explained by the refit (the sum of the mutations after refit minus the raw number of detected mutations before refit)
total_refit_variants <-  contri_tb2 %>% dplyr::group_by(Sample) %>% 
  dplyr::summarise(sum = sum(mean)) %>% 
  dplyr::ungroup()
total_refit_variants2 <- merge(total_refit_variants, Ageline_SBS_Counts_FA, by = "Sample")
total_refit_variants2$Unexplained <- total_refit_variants2$SBS_RAW - total_refit_variants2$sum
total_refit_variants2$Unexplained[total_refit_variants2$Unexplained < 0] <- NA
contri_tb3 <- rbind(contri_tb2, total_refit_variants2[,c("Sample", "Unexplained")] %>% as.data.frame() %>% 
  tidyr::gather(key = "sig", value = "mean", -Sample))


# Normalize the mutation counts for CallableLoci
contri_tb4 <- merge(contri_tb3, Metadata_Fig3[,c("Individual",  "Sample", "Age", "Phenotype", "CallableLoci_autosomal")], by = "Sample")
contri_tb4$mean_normalized <- contri_tb4$mean / as.numeric(contri_tb4$CallableLoci_autosomal)
contri_tb4 <- contri_tb4[which(contri_tb4$Age > 6 & contri_tb4$Age <20),]
contri_tb4 <- merge(contri_tb4, Ageline_SBS_Counts_FA[,c("Sample", "fit")], by = "Sample")

# Reorder the samples based on age and phenotype
contri_tb4$Sample <- factor(contri_tb4$Sample, levels = unique(contri_tb4$Sample[order(contri_tb4$Age)]))
contri_tb4$Individual <- factor(contri_tb4$Individual, levels = unique(contri_tb4$Individual[order(contri_tb4$Age)]))
contri_tb4$sig <- factor(contri_tb4$sig, levels = c("Unexplained", "PTA", "SBS3", "SBS18", "SBS1", "SBS5", "HSPC"))
contri_tb4$Phenotype <- factor(contri_tb4$Phenotype, levels = c("FanconiAnemia", "AML","Healthy"))

# For the healthy HSPCs we only want to show the mean contributions per individual (not per sample). Calculate the mean contributions per individual
means_per_individual <- contri_tb4 %>% dplyr::group_by(Individual, sig) %>% 
  dplyr::summarise(mean_normalized = mean(mean_normalized)) %>% 
  dplyr::ungroup()
means_per_individual$Sample <- means_per_individual$Individual
means_per_individual <- merge(means_per_individual, 
                              Metadata_Fig3[!duplicated(Metadata_Fig3$Individual),c("Individual", "Age", "Phenotype", "CallableLoci_autosomal")], 
                              by = "Individual", all.x = T)

contri_tb5 <- contri_tb4[contri_tb4$Phenotype != "Healthy",]
contri_tb5 <- rbind(contri_tb5[,c("Individual" ,"sig","mean_normalized","Sample","Age","Phenotype")], 
                    means_per_individual[means_per_individual$Phenotype == "Healthy",c("Individual" ,"sig","mean_normalized","Sample","Age","Phenotype")])
contri_tb5$sig <- factor(contri_tb5$sig, levels = c("Unexplained", "PTA", "SBS3", "SBS18", "SBS1", "SBS5", "HSPC"))
contri_tb5$Phenotype <- factor(contri_tb5$Phenotype, levels = c("FanconiAnemia", "AML","Healthy"))

Signature_colors <- c("#b3b3b3", "#fc8d62", "#e78ac3","#ffd92f", "#8da0cb","#66c2a5", "#a6d854"  )
names(Signature_colors) <- levels(contri_tb5$sig )

# A line will be added showing the expected number of SBSs based on age for each individual
contri_tb6 <- contri_tb4[!duplicated(contri_tb4$Sample),c("Sample", "Individual", "Phenotype", "fit")]
colnames(contri_tb6)[colnames(contri_tb6) == "fit"] <- "mean_normalized"

# An asterix will be added above the bars of the samples that are not amplified with PTA
# First calculate the total number of mutations per sample, to determine the y-position of the asterix
summed_mutations <- contri_tb5 %>% dplyr::group_by(Sample, Individual, Phenotype) %>% 
  dplyr::summarise(mean_normalized = sum(mean_normalized, na.rm = T))
summed_mutations$Label <- "*"
summed_mutations$Label[summed_mutations$Sample %in% Metadata_Fig3$Sample[Metadata_Fig3$Type == "PTA"]] <- ""

# The facet strips of the individuals will contain vertical text
facet_strips <- strip_nested(
  # Horizontal strips
  text_x = elem_list_text(angle = c(0,90), vjust = c(1,0.5), hjust = c(0.5,1)),
  by_layer_x = TRUE,
  )

ggplot(contri_tb5, aes(x = Sample, y = mean_normalized, fill = sig)) +
  geom_bar(stat = "identity") + 
  geom_text(data = summed_mutations, aes(x =Sample, label = Label, y = mean_normalized+20), inherit.aes = F, size = 1.6) +
  geom_hline(data = contri_tb6, aes(yintercept = mean_normalized), linewidth = 0.3) +
  facet_nested(.~Phenotype+Individual, scales = "free_x", space = "free_x", nest_line = element_line(),switch = "x", strip = facet_strips) +
  scale_fill_manual(values = Signature_colors) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1100)) +
  theme_classic(base_size = 6) +
  theme(axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        ggh4x.facet.nestline = element_line(colour = "black"),
        legend.key.height= unit(0.3, 'cm'),
        legend.key.width= unit(0.3, 'cm'),panel.spacing = unit(0.1, "lines"), legend.position = "bottom"
        ) +
  labs(fill = "Signature", y= "Base substitions" )

# There is still a large gap between the facet strips and the x-axis label (seems all strips, both vertical and horizontal have the same height) > adjust in illustrator
ggsave(paste(Output_dir_fig3, "3e_Fanconi_SignatureRefit_Abs.pdf", sep = ""), width = 6.7, height = 7, units = "cm",dpi = 300)

# Cosine similarities


mut_mat_individuals_rel

mut_mat_comparison <- cbind(mut_mat_grouped[,"Healthy"], signatures[,c("SBS1", "SBS3", "SBS5","SBS18", "PTA")])
colnames(mut_mat_comparison)[1] <- "Healthy\nHSCs"
# Select the mutation matrices from the FA patients
cos_sim <- melt(cos_sim_matrix(mut_mat_individuals_rel[,c(1:6)], mut_mat_comparison))

colnames(cos_sim) <- c("Individual", "MutProfile", "CosineSimilarity")
cos_sim$Individual <- factor(cos_sim$Individual, levels = unique(Metadata_FA$Individual[order(Metadata_FA$Age, decreasing = T)]))
ggplot(cos_sim, aes(x = MutProfile, y = Individual, fill = CosineSimilarity)) + geom_tile() +
  geom_text(aes(label = round(CosineSimilarity, 2)), size = 1.6) + 
  scale_fill_gradient(
    low = "cornsilk",
    high = "#D55E00",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill") +
  labs(fill = "Cosine\nSimilarity", x = "Mutational profiles", y = "Individual") +
  theme_classic(base_size = 6) +
  theme(        strip.background = element_blank(),
                strip.text.x = element_blank(),
                legend.key.height= unit(0.3, 'cm'),
                legend.key.width= unit(0.3, 'cm'),
                legend.position = "bottom",
                axis.text.x = element_text(angle = 60, hjust=1))

ggsave(paste(Output_dir_fig3, "3f_Fanconi_CosineSimilaritiesSignatures.pdf", sep = ""), width = 4.5, height = 5.5, units = "cm",dpi = 300)





### Supplemental data

Variant_Counts_m <- melt(Variant_Counts[,c("Sample", "SBS_RAW", "SBS_FAIL", "SBS_NORM", "Label2", "Individual", "Age", "fit", "Type")],
                         id.vars = c("Sample", "Label2", "Individual", "Age", "fit", "Type"))
Variant_Counts_m$Extrapolation <- ifelse(Variant_Counts_m$variable == "SBS_NORM", "E", "U")
Variant_Counts_m$Extrapolation <- factor(Variant_Counts_m$Extrapolation , levels = c("U", "E"))
Variant_Counts_m$Filter <- ifelse(Variant_Counts_m$variable == "SBS_FAIL", "FAIL", "PASS")
Variant_Counts_m$Individual <- factor(Variant_Counts_m$Individual, levels = unique(Variant_Counts_m$Individual[order(Variant_Counts_m$Age, decreasing = F)]))
Variant_Counts_m$Label2 <- ifelse(Variant_Counts_m$Type != "PTA", paste(Variant_Counts_m$Label2, "*", sep =""), Variant_Counts_m$Label2)
Variant_Counts_m_uncorrected <- Variant_Counts_m[Variant_Counts_m$variable != "SBS_NORM",]
Variant_Counts_uncorrected_freq <- Variant_Counts_m_uncorrected %>% dplyr::group_by(Individual, Label2, Filter) %>% 
  dplyr::summarise(sum = sum(value)) %>% 
  mutate(Freq = sum/sum(sum))


ggplot(Variant_Counts_m, aes( x = Extrapolation, y = value, fill = Filter)) + 
  geom_bar(stat = "identity", col = "white", linewidth = 0.2) + 
  geom_hline(aes(yintercept = fit)) +
  facet_nested(.~Individual + Label2, switch = "x") +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1500)) +
  theme_classic(base_size = 6) +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5), 
        strip.placement = "outside",
        legend.key.size = unit(0.25, 'cm'),
        legend.margin = margin(l = 0, r = 0),
        strip.background = element_rect(colour =  "white", fill = "gray95"),
        panel.spacing = unit(0.1, "lines")) +
  labs(y = "Single base substitions (#)", x = "Sample")

Output_dir_FigureS9 <- paste(Output_dir, "FigureS9/", sep = "")
if(dir.exists(Output_dir_FigureS9) == F){
  print(paste("# Creating dir: ", Output_dir_FigureS9, sep = ""))
  dir.create(Output_dir_FigureS9, recursive = T)
}
ggsave(paste(Output_dir_FigureS9, "/FigureS9A_Fanconi_PTATO_SBSfilter.pdf", sep = ""), width = 11, height = 5, units = "cm",dpi = 300)


ggplot(Variant_Counts_uncorrected_freq, aes( x = Label2, y = Freq, fill = Filter)) + 
  geom_bar(stat = "identity", col = "white", linewidth = 0.2) + 
  #geom_hline(aes(yintercept = fit)) +
  facet_nested(.~Individual , switch = "x", scales = "free", space = "free") +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1), labels = scales::percent) +
  theme_classic(base_size = 6) +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5), 
        strip.placement = "outside",
        legend.key.size = unit(0.25, 'cm'),
        strip.background = element_rect(colour =  "white", fill = "gray95"),
        strip.text = element_text(angle = 90),
        panel.spacing = unit(0.1, "lines")) +
  labs(y = "Single base substitions (%)", x = "Sample")

ggsave(paste(Output_dir_FigureS9, "/FigureS9B_Fanconi_PTATO_SBSfilter_rel.pdf", sep = ""), width = 8, height = 5, units = "cm",dpi = 300)
