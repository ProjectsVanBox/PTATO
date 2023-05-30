### Figure 5

library(RColorBrewer)
library(StructuralVariantAnnotation)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)
library(ggpubr)
library(dplyr)
library(rstatix)


options(scipen = 999)

input_dir <- "/path/to/MendeleyData/directory/"
output_dir_fig5 <- "/path/to/output/directory/"
if(dir.exists(output_dir_fig5) == F){dir.create(output_dir_fig5)}
output_dir_figS10 <- "/path/to/output/directory/"
if(dir.exists(output_dir_figS10) == F){dir.create(output_dir_figS10)}


# The raw VCFs contain the GRIPSS vcf output
SV_VCFs_RAW <- list("BULK" = paste(input_dir, "IBFM35/PTATO/intermediate/svs/gripss/IBFM35/IBFM35-DX2BM-MSCBULK/IBFM35-DX2BM-AMLBULK.gripss.somatic.filtered.vcf.gz", sep = ""),
                   PTA1 = paste(input_dir, "IBFM35/PTATO/intermediate/svs/gripss/IBFM35/IBFM35-DX2BM-MSCBULK/IBFM35-DX2BM-HSCPTAP1D9.gripss.somatic.filtered.vcf.gz", sep = ""), 
                   PTA2 = paste(input_dir, "IBFM35/PTATO/intermediate/svs/gripss/IBFM35/IBFM35-DX2BM-MSCBULK/IBFM35-DX2BM-HSCPTAP1E9.gripss.somatic.filtered.vcf.gz", sep = ""),
                   PTA3 = paste(input_dir, "IBFM35/PTATO/intermediate/svs/gripss/IBFM35/IBFM35-DX2BM-MSCBULK/IBFM35-DX2BM-HSCPTAP1G9.gripss.somatic.filtered.vcf.gz", sep = "")) 

# These are the PTATO-SV filtered VCFs:
SV_VCFs_Filtered <- list("BULK" = paste(input_dir, "IBFM35/SV/IBFM35-DX2BM-AMLBULK/purple/IBFM35-DX2BM-AMLBULKT.purple.sv.vcf.gz",  sep = ""),
                    PTA1 = paste(input_dir, "IBFM35/PTATO/intermediate/svs/Integration/IBFM35/IBFM35-DX2BM-MSCBULK/IBFM35-DX2BM-HSCPTAP1D9.integrated.svs.filtered.vcf",  sep = ""),
                    PTA2 = paste(input_dir, "IBFM35/PTATO/intermediate/svs/Integration/IBFM35/IBFM35-DX2BM-MSCBULK/IBFM35-DX2BM-HSCPTAP1E9.integrated.svs.filtered.vcf", sep = ""),
                    PTA3 = paste(input_dir, "IBFM35/PTATO/intermediate/svs/Integration/IBFM35/IBFM35-DX2BM-MSCBULK/IBFM35-DX2BM-HSCPTAP1G9.integrated.svs.filtered.vcf", sep = "")) 

# PTATO CNVs
PTATO_CNV_files <- list(PTA1 = paste(input_dir, "IBFM35/PTATO/intermediate/svs/Integration/IBFM35/IBFM35-DX2BM-MSCBULK/IBFM35-DX2BM-HSCPTAP1D9.integrated.cnvs.txt", sep = ""),
                        PTA2 = paste(input_dir, "IBFM35/PTATO/intermediate/svs/Integration/IBFM35/IBFM35-DX2BM-MSCBULK/IBFM35-DX2BM-HSCPTAP1E9.integrated.cnvs.txt", sep = ""),
                        PTA3 = paste(input_dir, "IBFM35/PTATO/intermediate/svs/Integration/IBFM35/IBFM35-DX2BM-MSCBULK/IBFM35-DX2BM-HSCPTAP1G9.integrated.cnvs.txt", sep = ""))
brewer.pal(n=10,"Paired")
SV_Colors <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FF7F00", "#CAB2D6","#6A3D9A")
names(SV_Colors) <- c("DUP <10kb", "DUP >10kb", "INS <10kb", "INS >10kb", "DEL <10kb", "DEL >10kb", "CTX", "INV <10kb","INV >10kb")

# Read and collect the GRIPSS and PTATO SV data
SV_Data <- data.frame()
for(Sample in names(SV_VCFs_RAW)){
  print(Sample)
  GRIPSS_file <- SV_VCFs_RAW[[Sample]]
  PTATO_file <- SV_VCFs_Filtered[[Sample]]
  
  print(GRIPSS_file)
  print(PTATO_file)
  
  GRIPSS <- readVcf(GRIPSS_file, "hg38")
  GRIPSS_PASS <- GRIPSS[rowRanges(GRIPSS)$FILTER == "PASS",]
  GRIPSS_GR <- breakpointRanges(GRIPSS_PASS)
  GRIPSS_PASS <- GRIPSS_PASS[names(GRIPSS_GR),]
  info(GRIPSS_PASS)$SVTYPE <- simpleEventType(GRIPSS_GR[names(GRIPSS_PASS)])
  info(GRIPSS_PASS)$SVLEN <- GRIPSS_GR[names(GRIPSS_PASS)]$svLen
  
  GRIPSS_PASS_f <- info(GRIPSS_PASS)$EVENT[-which(as.vector(seqnames(GRIPSS_PASS)) %in% c(1:22, "X", "Y"))]
  
  PTATO <- readVcf(PTATO_file, "hg38")
  PTATO_PASS <- PTATO[rowRanges(PTATO)$FILTER == "PASS",]
  PTATO_GR <- breakpointRanges(PTATO_PASS)
  PTATO_PASS <- PTATO_PASS[names(PTATO_GR),]
  info(PTATO_PASS)$SVTYPE <- simpleEventType(PTATO_GR[names(PTATO_PASS)])
  info(PTATO_PASS)$SVLEN <- PTATO_GR[names(PTATO_PASS)]$svLen
  
  GRIPSS_output <- info(GRIPSS_PASS)
  GRIPSS_output$Category <- ifelse(GRIPSS_output$SVLEN < 10000, "<10kb", ">10kb")
  GRIPSS_output$Label <- GRIPSS_output$SVTYPE
  GRIPSS_output$Data <- "Raw"
  GRIPSS_output$SV <- row.names(GRIPSS_output)
  # Count each SV event only once
  if(length(which(GRIPSS_output$EVENT %in% GRIPSS_PASS_f)) > 0){
    GRIPSS_output <- GRIPSS_output[-which(GRIPSS_output$EVENT %in% GRIPSS_PASS_f),]
  }
  GRIPSS_output <- GRIPSS_output[!duplicated(GRIPSS_output$EVENT),]
  
  PTATO_output <- info(PTATO_PASS)
  PTATO_output$Category <- ifelse(PTATO_output$SVLEN < 10000, "<10kb", ">10kb")
  PTATO_output$Label <- PTATO_output$SVTYPE
  PTATO_output$Data <- "Filtered"
  PTATO_output$SV <- row.names(PTATO_output)
  # Count each event only once
  PTATO_output <- PTATO_output[!duplicated(PTATO_output$EVENT),]
  
  output <- rbind(GRIPSS_output[,c("EVENT", "SV", "SVTYPE","SVLEN","Category","Label","Data")], PTATO_output[,c("EVENT", "SV", "SVTYPE","SVLEN","Category","Label","Data")])
  
  output$Label[output$SVTYPE == "DEL"] <- paste(output$Label[output$SVTYPE == "DEL"], output$Category[output$SVTYPE == "DEL"]  )
  output$Label[output$SVTYPE == "DUP"] <- paste(output$Label[output$SVTYPE == "DUP"], output$Category[output$SVTYPE == "DUP"]  )
  output$Label[output$SVTYPE == "INV"] <- paste(output$Label[output$SVTYPE == "INV"], output$Category[output$SVTYPE == "INV"]  )
  output$Label[output$SVTYPE == "INS"] <- paste(output$Label[output$SVTYPE == "INS"], output$Category[output$SVTYPE == "INS"]  )
  
  output$Sample <- Sample
  output <- as.data.frame(output)
  SV_Data <- rbind(SV_Data, output)
}

# Figure 5B (number of SV events)
SV_Data$Data <- factor(SV_Data$Data, levels = c("Raw", "Filtered"))
SV_Data$Data2 <- as.vector(SV_Data$Data )
SV_Data$Data2[SV_Data$Data2 == "Raw"] <- "Before PTATO"
SV_Data$Data2[SV_Data$Data2 == "Filtered"] <- "After PTATO"
SV_Data$Data2 <- factor(SV_Data$Data2, levels = c("Before PTATO", "After PTATO"))

ggplot(SV_Data, aes(x = Sample)) + 
  geom_bar(stat = "count", aes(fill = Label), width = 0.8) + 
  geom_text(stat='count', aes(label=after_stat(count), y =after_stat(count)), vjust=-0.5, size = 1.8) +
  facet_wrap(.~Data2) +
  scale_fill_manual(values = SV_Colors) + 
  theme_classic(base_size = 6) + 
  scale_y_continuous(expand = c(0,0.02), limits = c(0,650)) + 
  labs(y = "SV Events (Breakends)", fill = "SV Type") +
  theme(legend.key.size = unit(0.3, 'cm'), legend.margin = margin(0,0,0,0), panel.border = element_rect(color = "black", fill = NA), strip.background = element_rect(fill = "gray90"))

ggsave(paste(output_dir_fig5, "Figure5B.pdf", sep = ""), width = 6.8, height = 4, units = "cm")

# The breakpointRanges function does not work for the filtered SV VCf somehow. This is a work-around
PTA <- readVcf(paste(input_dir, "IBFM35/PTATO/intermediate/svs/Integration/IBFM35/IBFM35-DX2BM-MSCBULK/IBFM35-DX2BM-HSCPTAP1E9.integrated.svs.filtered.vcf", sep = ""), "hg38")
PTA_raw <- readVcf(paste(input_dir, "IBFM35/PTATO/intermediate/svs/gripss/IBFM35/IBFM35-DX2BM-MSCBULK/IBFM35-DX2BM-HSCPTAP1E9.gripss.somatic.filtered.vcf.gz", sep = ""), "hg38")

PTA_filtered <- PTA_raw[names(PTA_raw) %in% names(PTA),]
SVGR_PTA  <- breakpointRanges(PTA_filtered) 

Size_limit <- 10000

# Bulk (PURPLE output)
breakends_bulk_file <- paste(input_dir,"/IBFM35/SV/IBFM35-DX2BM-AMLBULK/purple/IBFM35-DX2BM-AMLBULKT.purple.sv.vcf.gz", sep = "")
breakends_bulk <- readVcf(breakends_bulk_file, "hg38")
breakends_bulk_filter <- breakends_bulk[rowRanges(breakends_bulk)$FILTER == "PASS"]
SVGR_bulk <- breakpointRanges(breakends_bulk_filter)
SVGR_bulk <- SVGR_bulk[which(abs(SVGR_bulk$svLen) > Size_limit | is.na(SVGR_bulk$svLen)),]

SVoverlap <- data.frame()

for(Sample in names(SV_VCFs_Filtered)[2:4]){
  print(Sample)
  PTATO_file <- SV_VCFs_Filtered[[Sample]]
  PTATO <- readVcf(PTATO_file, "hg38")
  PTATO_PASS <- PTATO[rowRanges(PTATO)$FILTER == "PASS",]
  PTATO_GR <- breakpointRanges(PTATO_PASS)
  PTATO_GR <- PTATO_GR[which(abs(PTATO_GR$svLen) > Size_limit | is.na(PTATO_GR$svLen)),]
  PTATO_GR <- PTATO_GR[as.character(seqnames(PTATO_GR)) %in% c(1:22),]
  #PTATO_PASS <- PTATO_PASS[names(PTATO_GR),]
  #info(PTATO_PASS)$SVTYPE <- simpleEventType(PTATO_GR[names(PTATO_PASS)])
  #info(PTATO_PASS)$SVLEN <- PTATO_GR[names(PTATO_PASS)]$svLen
  olap <- findBreakpointOverlaps(SVGR_bulk, PTATO_GR, maxgap = 1e5)
  
  overlapping <- SVGR_bulk[queryHits(olap)]
  print(overlapping)
  missing <- SVGR_bulk[-queryHits(olap),]
  SVoverlap_sample <- data.frame(Sample = Sample, Overlapping = length(unique(overlapping$event)), 
                                Missing = length(unique(missing$event)))
  SVoverlap_sample$Additional <- length(unique(PTATO_GR$event)) - SVoverlap_sample$Overlapping 
  SVoverlap <- rbind(SVoverlap, SVoverlap_sample)

  
}
SVoverlap_m <- melt(SVoverlap)

ggplot(SVoverlap_m, aes(x = Sample, y = value, fill = variable)) + geom_bar(stat = "identity")


# Collect purple CNVs for the bulk sample
PURPLE_CNV <- read.delim(paste(input_dir, "/IBFM35/SV/IBFM35-DX2BM-AMLBULK/purple/IBFM35-DX2BM-AMLBULKT.purple.cnv.somatic.tsv", sep = ""))
PURPLE_CNV2 <- PURPLE_CNV[PURPLE_CNV$copyNumber > 2.5 | PURPLE_CNV$copyNumber < 1.5,]

gain_purp <- GRanges(seqnames = PURPLE_CNV2$chromosome[PURPLE_CNV2$copyNumber > 2.5], IRanges(PURPLE_CNV2$start[PURPLE_CNV2$copyNumber > 2.5], PURPLE_CNV2$end[PURPLE_CNV2$copyNumber > 2.5]))
gain_purp2 <- reduce(gain_purp)
gain_purp3 <- gain_purp2[width(gain_purp2) > Size_limit,]
gain_purp3 <- gain_purp3[seqnames(gain_purp3) %in% c(1:22),]

loss_purp <- GRanges(seqnames = PURPLE_CNV2$chromosome[PURPLE_CNV2$copyNumber < 1.5], IRanges(PURPLE_CNV2$start[PURPLE_CNV2$copyNumber < 1.5], PURPLE_CNV2$end[PURPLE_CNV2$copyNumber < 1.5]))
loss_purp2 <- reduce(loss_purp)
loss_purp3 <- loss_purp2[width(loss_purp2) > Size_limit,]
loss_purp3 <- loss_purp3[seqnames(loss_purp3) %in% c(1:22),]



CNVoverlap <- data.frame()

for(Sample in names(PTATO_CNV_files)){
  print(Sample)
  PTATO_CNV <- read.delim(PTATO_CNV_files[[Sample]])
  PTATO_CNV <- PTATO_CNV[PTATO_CNV$chrom %in% c(1:22),]
  print(PTATO_CNV)
  
  PTATO_gains <- GRanges(PTATO_CNV$chrom[PTATO_CNV$CopyNumber == "Gain"], IRanges(PTATO_CNV$start[PTATO_CNV$CopyNumber == "Gain"],
                                                                                  PTATO_CNV$end[PTATO_CNV$CopyNumber == "Gain"]))
  olap_gains <- findOverlaps(gain_purp3, PTATO_gains, maxgap = 1e5)
  PTATO_losses <- GRanges(PTATO_CNV$chrom[PTATO_CNV$CopyNumber == "Loss"], IRanges(PTATO_CNV$start[PTATO_CNV$CopyNumber == "Loss"],
                                                                                  PTATO_CNV$end[PTATO_CNV$CopyNumber == "Loss"]))
  olap_losses <- findOverlaps(loss_purp3, PTATO_losses, maxgap = 1e5)
  
  
  CNVoverlap_sample <- data.frame(Sample = Sample, Overlapping = length(olap_gains) + length(olap_losses), 
                                 Missing = length(gain_purp3)+length(loss_purp3) - (length(olap_gains) + length(olap_losses)),
                                 Additional = length(PTATO_gains)+length(PTATO_losses) - (length(olap_gains) + length(olap_losses)))
  CNVoverlap <- rbind(CNVoverlap, CNVoverlap_sample)
}

# Merge SVs and CNVs
SVoverlap_m <- as.matrix(SVoverlap[,2:4])
row.names(SVoverlap_m) <- SVoverlap[,1]
CNVoverlap_m <- as.matrix(CNVoverlap[,2:4])
row.names(CNVoverlap_m) <- CNVoverlap[,1]

SV_counts <- SVoverlap_m + CNVoverlap_m
SV_counts_m <- melt(SV_counts)
SV_counts_m$Var2 <- factor(SV_counts_m$Var2, levels = c("Additional", "Missing", "Overlapping"))
ggplot(SV_counts_m, aes(x = Var1, y = value, fill = Var2)) + 
  geom_bar(stat = "identity") + 
  geom_text(aes(label = value), position = position_stack(vjust = 0.5), size = 2)+
  scale_y_continuous(expand = c(0,0)) +
  theme_classic(base_size = 6) +
  labs(y = "Structural variants (>10kb)", x= "Sample", fill = "Presence\nin bulk")+
  theme(legend.key.size = unit(0.4, 'cm'))

ggsave(paste(output_dir, "Fig5F_SharedSVs.png", sep = ""), width = 5, height = 4, dpi = 300, units = "cm")
ggsave(paste(output_dir, "Fig5F_SharedSVs.pdf", sep = ""), width = 5, height = 4, units = "cm")


### Figures 4D and 4E
# First bin the 1kb COBALT data of the bulk sample into 100kb windows for fair comparison with the PTA samples
IBFM35DX2BMAMLBULK_readcounts1kb <- read.delim(paste(input_dir, "/IBFM35/SV/IBFM35-DX2BM-AMLBULK/cobalt/IBFM35-DX2BM-AMLBULKT.cobalt.ratio.tsv", sep = ""))
IBFM35DX2BMAMLBULK_readcounts1kb$CopyNumber <- IBFM35DX2BMAMLBULK_readcounts1kb$tumorReadCount / mean(IBFM35DX2BMAMLBULK_readcounts1kb$tumorReadCount) 
CENTROMERE_FILE	<- paste(input_dir,"Resources/hg38_centromeres.txt", sep = "")
CENTROMERES <- read.delim(CENTROMERE_FILE)
CENTROMERES$chrom <- gsub("chr", "", CENTROMERES$chrom)
CENTROMERES_g <- GRanges(CENTROMERES$chrom, IRanges(CENTROMERES$chromStart, CENTROMERES$chromEnd))

RemoveCobaltOutliers <- function(COBALT_input, CENTROMERES = NULL, COBALT_MIN_COV = 0.01, COBALT_MAX_COV = 0.99){
  # Filter the bins with extreme high or low PON_meanReadCount
  COBALT_input$FILTER <- ""
  COBALT_input$FILTER[COBALT_input$PON_meanReadCount < quantile(COBALT_input$PON_meanReadCount, COBALT_MIN_COV)] <- "LOW_MEAN_COUNT"
  COBALT_input$FILTER[COBALT_input$PON_meanReadCount > quantile(COBALT_input$PON_meanReadCount, COBALT_MAX_COV)] <- "HIGH_MEAN_COUNT"
  
  if(!is.null(CENTROMERES)){
    if(class(CENTROMERES) == "GRanges"){
      # Overlap the COBALT bins with the centromeres
      COBALT_g <- GRanges(seqnames = COBALT_input$chromosome, IRanges(start = COBALT_input$position, 
                                                                      end = COBALT_input$position + 999))
      olap_COBALT_CENT <- findOverlaps(COBALT_g, CENTROMERES, maxgap = 5e6)
      
      COBALT_input$FILTER[queryHits(olap_COBALT_CENT)][COBALT_input$FILTER[queryHits(olap_COBALT_CENT)] != ""] <- paste(COBALT_input$FILTER[queryHits(olap_COBALT_CENT)][COBALT_input$FILTER[queryHits(olap_COBALT_CENT)] != ""], "CENTROMERIC", sep = ";")
      COBALT_input$FILTER[queryHits(olap_COBALT_CENT)][COBALT_input$FILTER[queryHits(olap_COBALT_CENT)] == ""] <- "CENTROMERIC"
    }
  }
  
  COBALT_input$FILTER[COBALT_input$FILTER == ""] <- "PASS"
  return(COBALT_input)
}

IBFM35DX2BMAMLBULK_readcounts1kb_ann <- RemoveCobaltOutliers(IBFM35DX2BMAMLBULK_readcounts1kb, CENTROMERES = CENTROMERES_g)
IBFM35DX2BMAMLBULK_readcounts1kb_filtered <- IBFM35DX2BMAMLBULK_readcounts1kb_ann[IBFM35DX2BMAMLBULK_readcounts1kb_ann$FILTER == "PASS",]

MergeBins <- function(Input,  Value_column = 3, Binsize = 100000, Ploidy = 2) {
  df <- data.frame()
  # Use mixedsort to order the chromosomes alphanumerically (allows different ref genomes, but causes dependency on gtools)
  Chromosomes <- gtools::mixedsort(as.character(unique((Input[,1]))))
  #print(Chromosomes)
  for (i in c(1:length(Chromosomes))) {
    Chrom <- Chromosomes[i]
    #print(Chrom)
    # Select chromosome
    tmp <- subset( Input, chromosome==Chrom)
    
    # Binning of positions
    grouping <- cut(tmp$position, 
                    c(seq(1,max(tmp$position),Binsize), max(tmp$position)), 
                    labels=seq(1,max(tmp$position),Binsize), max(tmp$position))
    
    # Calculate median readdepth
    medians <- data.frame(tapply(tmp[,Value_column], grouping, median))
    rm(tmp)
    
    colnames(medians) <- "medianReadCount"
    
    # Create clean data farme
    medians$Chromosome <- Chrom
    medians$start.pos <- as.numeric(as.character(rownames(medians)))
    medians$end.pos <- as.numeric(as.character(rownames(medians)))+Binsize-1
    # # Limit end of bins to end of chromosome
    # if(!is.null(Chr_sizes)){
    #   medians$end.pos[medians$end.pos >  as.numeric(Chr_sizes[Chrom])] <-  as.numeric(Chr_sizes[Chrom])
    # }
    medians$CopyNumber <- medians$medianReadCount*Ploidy
    if(length(which(is.na(medians$medianReadCount))) > 0){
      medians <- medians[-which(is.na(medians$medianReadCount)),]
    }
    # Merge data
    df <- rbind(df, medians)
  }
  df <- df[,c("Chromosome","start.pos", "end.pos", "medianReadCount", "CopyNumber")]
  return(df)
}

IBFM35DX2BMAMLBULK_readcounts100kb <- MergeBins(IBFM35DX2BMAMLBULK_readcounts1kb_filtered , Binsize = 100000, Value_column = 8)
IBFM35DX2BMAMLBULK_readcounts100kb$Sample <- "BULK"
IBFM35DX2BMAMLBULK_readcounts100kb$CopyNumber[IBFM35DX2BMAMLBULK_readcounts100kb$CopyNumber > 4.5] <- 4.5

# Read the 100kb filtered COBALT files of the PTA samples
readCounts_files <- list(PTA1 = paste(input_dir, "IBFM35/PTATO/intermediate/cnvs/readCounts/IBFM35-DX2BM-HSCPTAP1D9/IBFM35-DX2BM-HSCPTAP1D9.readcounts.filtered.100kb.txt", sep = ""), 
                         PTA2 = paste(input_dir, "IBFM35/PTATO/intermediate/cnvs/readCounts/IBFM35-DX2BM-HSCPTAP1E9/IBFM35-DX2BM-HSCPTAP1E9.readcounts.filtered.100kb.txt", sep = ""),
                         PTA3 = paste(input_dir, "IBFM35/PTATO/intermediate/cnvs/readCounts/IBFM35-DX2BM-HSCPTAP1G9/IBFM35-DX2BM-HSCPTAP1G9.readcounts.filtered.100kb.txt", sep = ""))

readCounts_overview <- data.frame()
for(sample in names(readCounts_files)){
  print(sample)
  print(readCounts_files[[sample]])
  readCounts <- read.delim(readCounts_files[[sample]])
  readCounts <- readCounts[,-5]
  colnames(readCounts)[4] <- "CopyNumber"
  readCounts$CopyNumber[readCounts$CopyNumber > 4.5] <- 4.5
  readCounts$Sample <- sample
  readCounts_overview <- rbind(readCounts_overview, readCounts)
}
readCounts_overview <- rbind(readCounts_overview, IBFM35DX2BMAMLBULK_readcounts100kb[,c(1:3, 5,6)])
readCounts_overview$Chromosome <- factor(readCounts_overview$Chromosome, levels = c(1:22, "X", "Y"))

# Read the segment data
readCountsSegments_files <- list(PTA1 = paste(input_dir, "IBFM35/PTATO/intermediate/cnvs/readCounts/IBFM35-DX2BM-HSCPTAP1D9/IBFM35-DX2BM-HSCPTAP1D9.readcounts.segments.txt", sep = ""), 
                                 PTA2 = paste(input_dir, "IBFM35/PTATO/intermediate/cnvs/readCounts/IBFM35-DX2BM-HSCPTAP1E9/IBFM35-DX2BM-HSCPTAP1E9.readcounts.segments.txt", sep = ""),
                                 PTA3 = paste(input_dir, "IBFM35/PTATO/intermediate/cnvs/readCounts/IBFM35-DX2BM-HSCPTAP1G9/IBFM35-DX2BM-HSCPTAP1G9.readcounts.segments.txt", sep = ""))

readCountsSegments_overview <- data.frame()
for(sample in names(readCountsSegments_files)){
  print(sample)
  print(readCountsSegments_files[[sample]])
  readCountsSegments <- read.delim(readCountsSegments_files[[sample]])[,1:4]
  names(readCountsSegments) <- c("Chromosome", "start.pos", "end.pos", "CopyNumber")
  readCountsSegments$CopyNumber[readCountsSegments$CopyNumber > 4.5] <- 4.5
  readCountsSegments$Sample <- sample
  readCountsSegments_overview <- rbind(readCountsSegments_overview, readCountsSegments)
} 

# Segments bulk sample:
IBFM35DX2BMAMLBULK_readcounts_segments <- read.delim(paste(input_dir, "IBFM35/SV/IBFM35-DX2BM-AMLBULK/cobalt/IBFM35-DX2BM-AMLBULKT.cobalt.ratio.pcf", sep = ""))
names(IBFM35DX2BMAMLBULK_readcounts_segments) <- c( "sampleID", "Chromosome","arm", "start.pos","end.pos","n.probes", "mean")
IBFM35DX2BMAMLBULK_readcounts_segments$CopyNumber <- 2 + IBFM35DX2BMAMLBULK_readcounts_segments$mean *2
IBFM35DX2BMAMLBULK_readcounts_segments$CopyNumber[IBFM35DX2BMAMLBULK_readcounts_segments$mean < 0] <- 2 + IBFM35DX2BMAMLBULK_readcounts_segments$mean[IBFM35DX2BMAMLBULK_readcounts_segments$mean < 0]

IBFM35DX2BMAMLBULK_readcounts_segments$Sample <- "BULK"
IBFM35DX2BMAMLBULK_readcounts_segments$CopyNumber[IBFM35DX2BMAMLBULK_readcounts_segments$CopyNumber > 4.5] <- 4.5
readCountsSegments_overview <- rbind(readCountsSegments_overview, IBFM35DX2BMAMLBULK_readcounts_segments[,c("Chromosome", "start.pos","end.pos","CopyNumber","Sample")])

readCountsSegments_overview$Chromosome <- factor(readCountsSegments_overview$Chromosome, levels = c(1:22, "X", "Y"))

# CNV calls
# CNVs bulk
BULK_cnvs <- read.delim(paste(input_dir,"/IBFM35/SV/IBFM35-DX2BM-AMLBULK/purple/circos/IBFM35-DX2BM-AMLBULKT.cnv.circos", sep =""))
BULK_cnvs <- BULK_cnvs[BULK_cnvs$value < -0.8 | BULK_cnvs$value > 0.8,]
BULK_cnvs$CopyNumber <- ifelse(BULK_cnvs$value < -1, "Loss", "Gain")
BULK_cnvs$Sample <-  "BULK"
names(BULK_cnvs) <- c("Chromosome", "start.pos", "end.pos", "value", "CopyNumber", "Sample")
BULK_cnvs[,1] <- gsub(pattern = "hs", replacement = "", x = BULK_cnvs[,1])

PURPLE_CNV3 <- PURPLE_CNV[PURPLE_CNV$copyNumber < 1.5 | PURPLE_CNV$copyNumber > 2.5 | PURPLE_CNV$baf > 0.8,]
PURPLE_CNV3 <- PURPLE_CNV3[PURPLE_CNV3$bafCount > 100 | PURPLE_CNV3$copyNumber <1.5,]
BULK_cnvs2 <- PURPLE_CNV3[,1:3]
BULK_cnvs2$CopyNumber <- "cnLOH"
BULK_cnvs2$CopyNumber[PURPLE_CNV3$copyNumber < 1.5] <- "Loss"
BULK_cnvs2$CopyNumber[PURPLE_CNV3$copyNumber > 2.5] <- "Gain"
BULK_cnvs2$Sample <- "BULK"
names(BULK_cnvs2) <- c("Chromosome", "start.pos", "end.pos", "CopyNumber", "Sample")

CNV_files <- list(PTA1 = paste(input_dir, "IBFM35/PTATO/intermediate/svs/Integration/IBFM35/IBFM35-DX2BM-MSCBULK/IBFM35-DX2BM-HSCPTAP1D9.integrated.cnvs.txt", sep = ""), 
                  PTA2 = paste(input_dir, "IBFM35/PTATO/intermediate/svs/Integration/IBFM35/IBFM35-DX2BM-MSCBULK/IBFM35-DX2BM-HSCPTAP1E9.integrated.cnvs.txt", sep = ""),
                  PTA3 = paste(input_dir, "IBFM35/PTATO/intermediate/svs/Integration/IBFM35/IBFM35-DX2BM-MSCBULK/IBFM35-DX2BM-HSCPTAP1G9.integrated.cnvs.txt", sep = ""))

CNV_overview <- data.frame()

for(sample in names(CNV_files)){
  print(sample)
  print(CNV_files[[sample]])
  CNVs <- read.delim(CNV_files[[sample]])
  CNVs$Sample <- sample
  CNV_overview <- rbind(CNV_overview, CNVs[,c(1,2,3,6,which(names(CNVs) == "Sample"))])
}
names(CNV_overview) <- c("Chromosome", "start.pos", "end.pos", "CopyNumber", "Sample")
CNV_overview <- rbind(CNV_overview, BULK_cnvs2)
CNV_overview$Chromosome <- factor(CNV_overview$Chromosome, levels = c(1:22, "X", "Y"))
CNV_overview$CopyNumber[CNV_overview$CopyNumber == "LOH"] <- "cnLOH"
### Color settings
sample_colors <- c( "honeydew4", "darkgrey","honeydew4", "darkgrey")
names(sample_colors) <- names(CNV_files)

CNV_cols <- c("steelblue3","red3","darkorange", "grey60")
names(CNV_cols) <- c("Gain","Loss","cnLOH", "Normal")


### Figure 5D
p <- ggplot(data = readCounts_overview[readCounts_overview$Chromosome %in% c(1:22),]) +
  geom_rect(data = CNV_overview[CNV_overview$Chromosome %in% c(1:22),], mapping = aes(xmin = start.pos, xmax = end.pos, fill = CopyNumber), ymin = -Inf, ymax = Inf, alpha = 0.25) +
  geom_point(data = readCounts_overview[readCounts_overview$Chromosome %in% c(1:22),], aes(x = (start.pos+end.pos) / 2, y = CopyNumber, col = Sample), alpha = 0.3, size = 0.01, show.legend=FALSE) +
  geom_segment(data = readCountsSegments_overview[readCountsSegments_overview$Chromosome %in% c(1:22),], 
               aes(x = start.pos, xend = end.pos, y = CopyNumber, yend = CopyNumber), col = "black") +
  #geom_hline(yintercept = 0) +
  facet_grid(Sample~Chromosome, space = "free_x", scales = "free", switch="both") +
  theme_bw(base_size = 5) +
  guides(colour = guide_legend(override.aes = list(size=1))) +
  labs(x = "Chromosome") +
  scale_color_manual(values = sample_colors) +
  scale_fill_manual(values = CNV_cols) +
  scale_y_continuous(limits = c(0,4.5)) +
  theme(
    axis.line.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.grid.minor=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.major.y=element_line(colour="grey90", linetype=3, size=0.1),
    panel.spacing = unit(0.05, "lines"),
    plot.title = element_text(hjust = 0.5),
    legend.position="bottom",
    panel.background=element_blank(),
    panel.border =element_blank(),
    strip.placement = "outside", legend.key.size = unit(0.4, 'cm'))

p

q <- ggplotGrob(p)
lg <- linesGrob(x=unit(c(0,0),"npc"), y=unit(c(0,1),"npc"),
                gp=gpar(col="gray", lwd=1))

for (k in grep("strip-b",q$layout$name)) {
  q$grobs[[k]]$grobs[[1]]$children[[1]] <- lg
}

ggsave(plot = arrangeGrob(q), paste(output_dir_fig5, "Fig5D_IBFM35_CopyNumber.png", sep = ""), width = 9, height = 7, dpi = 600, units = "cm")
ggsave(plot = arrangeGrob(q), paste(output_dir_fig5, "Fig5D_IBFM35_CopyNumber.pdf", sep = ""), width = 9, height = 7, units = "cm")

### Figure 5E

BAF_files <- list("BULK" = paste(input_dir, "IBFM35/PTATO/intermediate/cnvs/BAF/IBFM35/IBFM35-DX2BM-MSCBULK/IBFM35-DX2BM-AMLBULK.baf.binned.100kb.txt", sep = ""),
                  PTA1 = paste(input_dir, "IBFM35/PTATO/intermediate/cnvs/BAF/IBFM35/IBFM35-DX2BM-MSCBULK/IBFM35-DX2BM-HSCPTAP1D9.baf.binned.100kb.txt", sep = ""), 
                  PTA2 = paste(input_dir, "IBFM35/PTATO/intermediate/cnvs/BAF/IBFM35/IBFM35-DX2BM-MSCBULK/IBFM35-DX2BM-HSCPTAP1E9.baf.binned.100kb.txt", sep = ""),
                  PTA3 = paste(input_dir, "IBFM35/PTATO/intermediate/cnvs/BAF/IBFM35/IBFM35-DX2BM-MSCBULK/IBFM35-DX2BM-HSCPTAP1G9.baf.binned.100kb.txt", sep = ""))

BAF_overview <- data.frame()

for(sample in names(BAF_files)){
  print(sample)
  print(BAF_files[[sample]])
  BAFs <- read.delim(BAF_files[[sample]])
  BAFs$Sample <- sample
  BAF_overview <- rbind(BAF_overview, BAFs)
}
BAF_overview$Chromosome <- factor(BAF_overview$Chromosome, levels = c(1:22, "X", "Y"))

#
BAFSegment_files <- list("BULK" = paste(input_dir, "IBFM35/PTATO/intermediate/cnvs/BAF/IBFM35/IBFM35-DX2BM-MSCBULK/IBFM35-DX2BM-AMLBULK.baf.segments.txt", sep = ""),
                         PTA1 = paste(input_dir, "IBFM35/PTATO/intermediate/cnvs/BAF/IBFM35/IBFM35-DX2BM-MSCBULK/IBFM35-DX2BM-HSCPTAP1D9.baf.segments.txt", sep = ""), 
                         PTA2 = paste(input_dir, "IBFM35/PTATO/intermediate/cnvs/BAF/IBFM35/IBFM35-DX2BM-MSCBULK/IBFM35-DX2BM-HSCPTAP1E9.baf.segments.txt", sep = ""),
                         PTA3 = paste(input_dir, "IBFM35/PTATO/intermediate/cnvs/BAF/IBFM35/IBFM35-DX2BM-MSCBULK/IBFM35-DX2BM-HSCPTAP1G9.baf.segments.txt", sep = ""))

BAFSegments_overview <- data.frame()

for(sample in names(BAFSegment_files)){
  print(sample)
  print(BAFSegment_files[[sample]])
  BAFSegments <- read.delim(BAFSegment_files[[sample]])
  names(BAFSegments) <- names(BAF_overview)[1:4]
  BAFSegments$Sample <- sample
  BAFSegments_overview <- rbind(BAFSegments_overview, BAFSegments)
}
BAFSegments_overview$Chromosome <- factor(BAFSegments_overview$Chromosome, levels = c(1:22, "X", "Y"))


p2 <- ggplot(data = readCounts_overview[readCounts_overview$Chromosome %in% c(1:22),]) +
  geom_rect(data = CNV_overview[CNV_overview$Chromosome %in% c(1:22),], mapping = aes(xmin = start.pos, xmax = end.pos, fill = CopyNumber), ymin = -Inf, ymax = Inf, alpha = 0.25) +
  geom_point(data = BAF_overview[BAF_overview$Chromosome %in% c(1:22),], aes(x = (start.pos+end.pos) / 2, y = medianBAF, col = Sample), alpha = 0.3, size = 0.01, show.legend=FALSE) +
  geom_segment(data = BAFSegments_overview[BAFSegments_overview$Chromosome %in% c(1:22),],  aes(x = start.pos, xend = end.pos, y = medianBAF, yend = medianBAF), col = "black") +
  #geom_hline(yintercept = 0) +
  facet_grid(Sample~Chromosome, space = "free_x", scales = "free", switch="both") +
  theme_bw(base_size = 5) +
  guides(colour = guide_legend(override.aes = list(size=1))) +
  scale_y_reverse() +
  labs(x = "Chromosome", y = "Deviation of allele frequency (DAF)") +
  scale_color_manual(values = sample_colors) +
  scale_fill_manual(values = CNV_cols) +
  theme(
    axis.line.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.grid.minor=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.major.y=element_line(colour="grey90", linetype=3, size=0.1),
    panel.spacing = unit(0.05, "lines"),
    plot.title = element_text(hjust = 0.5),
    legend.position="bottom",
    panel.background=element_blank(),
    panel.border =element_blank(),
    strip.placement = "outside", legend.key.size = unit(0.4, 'cm'))

p2

q2 <- ggplotGrob(p2)
lg <- linesGrob(x=unit(c(0,0),"npc"), y=unit(c(0,1),"npc"),
                gp=gpar(col="gray", lwd=1))

for (k in grep("strip-b",q$layout$name)) {
  q2$grobs[[k]]$grobs[[1]]$children[[1]] <- lg
}

ggsave(plot = arrangeGrob(q2), paste(output_dir_fig5, "Fig5E_IBFM35_DAF.png", sep = ""), width = 9, height = 7, dpi = 600, units = "cm")
ggsave(plot = arrangeGrob(q2), paste(output_dir_fig5, "Fig5E_IBFM35_DAF.pdf", sep = ""), width = 9, height = 7, units = "cm")


### Figures 5G-5H
SV_Counts <- read.delim(paste(input_dir, "SV_Counts_FA.txt", sep = ""))
SV_Counts <- SV_Counts[!is.na(SV_Counts$Large.deletions),]
SV_Counts$Type[SV_Counts$Type == "FanconiAnemia"] <- "Fanconi\nAnemia"
SV_Counts$Type[SV_Counts$Type == "Healthy_HSPC"] <- "Healthy"
SV_Counts$Type <- factor(SV_Counts$Type, levels = c("Healthy", "Fanconi\nAnemia"))

DeletionSummary <- as.data.frame(SV_Counts %>% group_by(Type) %>%
                                   summarise(across(Large.deletions, 
                                                    list(mean = mean, min = min, max = max, sd = sd, group_size = length))))

# Showing multiple testing adjusted p-values does not work for the standard stat_compare_means function
wilcox <- compare_means(Large.deletions ~ Type,  data = SV_Counts)
comparisons_with_healthy <- wilcox[wilcox$group1 == "Healthy",]
# Only show the p-values for the significant differences between healthy and patients
my_comparisons <- list( c("Healthy", "PMCFANC03"), c("Healthy", "PMCFANC02"), c("Healthy", "IBFM35") )

stat.test <- wilcox_test(Large.deletions ~ Type, data = SV_Counts) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "Type", dodge = 1, step.increase = .12)
stat.test$Type <- stat.test$group2
DeletionSummary$Type[DeletionSummary$Type == "FanconiAnemia"] <- "Fanconi\nAnemia"
DeletionSummary$Type[DeletionSummary$Type == "Healthy_HSPC"] <- "Healthy"
DeletionSummary$Type <- factor(DeletionSummary$Type, levels = c("Healthy", "Fanconi\nAnemia"))
stat.test$y.position <- stat.test$y.position + 0.8


ggplot(DeletionSummary, aes(x = Type, y = Large.deletions_mean)) +
  geom_bar(stat = "identity", col = "black", fill = "gray", width = 0.8, linewidth = 0.2) +
  geom_errorbar(aes(ymin=Large.deletions_mean, ymax=Large.deletions_mean+Large.deletions_sd), linewidth = 0.2, width = 0.2,
                position=position_dodge(0.05)) +
  geom_text(aes(x =Type , label = paste("n=", Large.deletions_group_size, sep = "")), y = 2.6, size = 1.8, fontface = "italic") +
  geom_jitter(data = SV_Counts, aes(x = Type, y = Large.deletions, color = Sequencing), height = 0, width = 0.3, size = 0.5, alpha = 0.8) +
  stat_pvalue_manual(
    stat.test, label = "P={p.adj}", size = 1.6, bracket.size = 0.2, vjust = -0.3, tip.length = 0.03) +
  scale_color_brewer(palette = "Dark2") +
  coord_cartesian(ylim = c(0,3)) +
  scale_y_continuous(limits = c(0,3), expand = c(0,0)) +
  theme_classic(base_size = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.25, 'cm')) +
  labs(y = "Large deletions\n(per genome)", color = "Method", x = "Type of HSPC")
ggsave(paste(output_dir_fig5, "Figure5G.pdf", sep =""), width = 4, height = 4, units = "cm")

SVs_FA <- read.delim(paste(input_dir, "Table_S3.txt", sep = ""))
Deletions_FA <- SVs_FA[SVs_FA$SVType == "DEL",]
Deletions_FA$Type[Deletions_FA$Type == "FanconiAnemia"] <- "Fanconi\nAnemia"
Deletions_FA$Type[Deletions_FA$Type == "Healthy_HSPC"] <- "Healthy"
Deletions_FA$Type <- factor(Deletions_FA$Type, levels = c("Healthy", "Fanconi\nAnemia"))

Deletions_FA_Summary <- as.data.frame(Deletions_FA %>% group_by(Type) %>%
                                        summarise(across(Size, 
                                                         list(mean = mean, median = median, min = min, max = max, sd = sd, group_size = length))))
stat.test <- wilcox_test(Size ~ Type, data = Deletions_FA) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "Type", dodge = 1, step.increase = .12)

ggplot(Deletions_FA, aes(x = Type, y = Size)) +
  geom_boxplot(outlier.shape = NA, fill = "gray", linewidth = 0.2) + 
  geom_jitter(data = Deletions_FA, aes(x = Type, y = Size, color = Sequencing), height = 0, width = 0.3, size = 0.5, alpha = 0.8) +
  geom_text(data = Deletions_FA_Summary, aes(x = Type, label = paste("n=", Size_group_size, sep = ""), y = 80000), size = 1.8, fontface = "italic", inherit.aes = F) +
  scale_y_log10(breaks = c(10,100,1000,10000,100000), labels = c(10,100,1000,10000,100000))  +
  scale_color_brewer(palette = "Dark2") +
  coord_cartesian(ylim = c(10,100000)) +
  theme_classic(base_size = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.25, 'cm')) +
  labs(y = "Deletion size\n(basepairs)", color = "Method", x = "Type of HSPC")
ggsave(paste(output_dir_fig5, "Figure5H_FA_DeletionSizes.pdf", sep =""), width = 4.8, height = 4, units = "cm")


### Figure 5I

indel_sizes <- data.frame()

for(individual in unique(Metadata_Fig5$Individual)){
  print(individual)
  for(sample in Metadata_Fig5$Sample[Metadata_Fig5$Individual == individual]){
    print(sample)
    indels_sample <- indels_autosomal[[sample]]
    indel_size_sample <- data.frame(Individual = individual, Sample = sample, 
                                    Type = ifelse(lengths(indels_sample$REF) > lengths(unlist(indels_sample$ALT)), "DEL", "INS"),
                                    size = abs(lengths(indels_sample$REF) - lengths(unlist(indels_sample$ALT))))
    if(nrow(indel_sizes) > 0 ){
      indel_sizes <- rbind(indel_sizes, indel_size_sample)
    } else {
      indel_sizes <- indel_size_sample
    }
  }
}

indel_sizes$Category <- ifelse(indel_sizes$size > 1, "2-10bp", "1bp")
indel_sizes$Category <- ifelse(indel_sizes$size > 10, "11-50bp", indel_sizes$Category)
indel_sizes$Category <- ifelse(indel_sizes$size > 50, ">50bp", indel_sizes$Category)
indel_sizes$Category <- factor(indel_sizes$Category, levels = c(">50bp", "11-50bp","2-10bp", "1bp"))

# PEr sample
indel_sizes_PerSample <- indel_sizes[indel_sizes$Type == "DEL",] %>% group_by(Sample, Category, .drop = F) %>% 
  summarise(Freq = n()) %>%
  mutate(Freq = Freq/sum(Freq))
indel_sizes_PerSample <- merge(indel_sizes_PerSample, Metadata[!duplicated(Metadata$Individual),c("Sample", "Individual", "Age", "Phenotype")], by = "Sample")

# mean per individual
indel_sizes_perIndividual <- indel_sizes_PerSample %>% 
  group_by(Individual, Category) %>% 
  summarise_at(vars(Freq), list(mean = mean))
indel_sizes_perIndividual <- merge(indel_sizes_perIndividual, 
                                   Metadata[!duplicated(Metadata$Individual),c("Individual", "Age", "Phenotype")], by = "Individual")

# Per phenotype

indel_sizes_perGroup <- as.data.frame(indel_sizes_perIndividual %>% 
                                        group_by(Phenotype, Category)  %>% 
                                        summarise_at(vars(mean), list(mean = mean)))
indel_sizes_perGroup$Species <- "Human"
indel_sizes_Mice <- data.frame(Phenotype = c(rep("Fancd2",4),
                                             rep("Aldh2;Fancd2", 4)), 
                               Category = rep(c(">50bp", "11-50bp", "2-10bp", "1bp"),2),
                               mean = c(0, 0.286, 0.171, 0.543,
                                        0.172, 0.303, 0.182, 0.343),
                               Species = "Mouse")

indel_sizes_FA <- rbind(indel_sizes_perGroup, indel_sizes_Mice)
indel_sizes_FA$Phenotype <- factor(indel_sizes_FA$Phenotype, levels = c("Healthy", "FanconiAnemia", "AML", "Fancd2", "Aldh2;Fancd2"))

ggplot(indel_sizes_perGroup, aes(x = Phenotype, y = mean)) + geom_bar(stat = "identity", aes(fill = Category), col = "black") +
  #geom_text(data = indel_sizes_perSample, aes(label = Freq, y = 1.05, x = Sample))+
  theme_classic(base_size  = 8) + scale_fill_brewer(palette = "Reds", direction = -1) +
  labs(y = "Relative contribution")



ggplot(indel_sizes_perIndividual, aes(x = Individual, y = mean)) + geom_bar(stat = "identity", aes(fill = Category), col = "black") +
  #geom_text(data = indel_sizes_perSample, aes(label = Freq, y = 1.05, x = Sample))+
  theme_classic(base_size  = 8) + scale_fill_brewer(palette = "Reds", direction = -1) +
  labs(y = "Relative contribution")

ggplot(indel_sizes_FA, aes(x = Phenotype, y = mean)) + 
  geom_bar(stat = "identity", aes(fill = Category), col = "black", size = 0.2) +
  facet_grid(.~Species, scales = "free", space = "free") +
  #geom_text(data = indel_sizes_perSample, aes(label = Freq, y = 1.05, x = Sample))+
  theme_classic(base_size  = 6) + 
  scale_fill_brewer(palette = "Reds", direction = -1) +
  scale_y_continuous(limits = c(0, 1.1), expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), 
        legend.key.size = unit(0.2, 'cm'),
        legend.position = "bottom") +
  labs(y = "Relative contribution", fill = "Deletion size")

ggsave(paste(output_dir_fig5, "5I_Fanconi_IndelSizes.pdf", sep = ""), width = 5, height = 5, units = "cm",dpi = 300)

