# Figure S11C and S11D: Copy number variants in AHH-1 cell lines and IBFM26

library(ggplot2)
library(GenomicRanges)
library(grid)
library(gridExtra)

options(scipen = 999)
Input_dir <- "/path/to/MendeleyData/directory/"
Output_dir <- "/path/to/output/directory/"

# First bin the 1kb COBALT data of the bulk-wgs sample into 100kb windows for fair comparison with the PTA samples
BULK_readcounts1kb <- read.delim(paste(Input_dir, "PMCAHH1-WT/SV/PMCAHH1-WT-C6/cobalt/PMCAHH1-WT-C6T.cobalt.ratio.tsv", sep = ""))
BULK_readcounts1kb$medianReadcount <- BULK_readcounts1kb$tumorReadCount / mean(BULK_readcounts1kb$tumorReadCount) * 2
BULK_readcounts1kb$CopyNumber <- BULK_readcounts1kb$medianReadcount 

CENTROMERE_FILE	<- paste(Input_dir, "Resources/hg38_centromeres.txt", sep = "")

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

BULK_readcounts1kb_ann <- RemoveCobaltOutliers(BULK_readcounts1kb, CENTROMERES = CENTROMERES_g)
BULK_readcounts1kb_filtered <- BULK_readcounts1kb_ann[BULK_readcounts1kb_ann$FILTER == "PASS",]

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

BULK_readcounts100kb <- MergeBins(BULK_readcounts1kb_filtered , Binsize = 100000, Value_column = 8)
BULK_readcounts100kb$Sample <- "CLONE_C6"
BULK_readcounts100kb$CopyNumber[BULK_readcounts100kb$CopyNumber > 4.5] <- 4.5

# Segments bulk
BULK_readcounts_segments <- read.delim(paste(Input_dir,"PMCAHH1-WT/SV/PMCAHH1-WT-C6/cobalt/PMCAHH1-WT-C6T.cobalt.ratio.pcf", sep = ""))
BULK_readcounts_segments$Sample <-  "CLONE_C6"
colnames(BULK_readcounts_segments)[colnames(BULK_readcounts_segments) == "mean"] <- "CopyNumber"
colnames(BULK_readcounts_segments)[colnames(BULK_readcounts_segments) == "chrom"] <- "Chromosome"
BULK_readcounts_segments$CopyNumber <- 2+BULK_readcounts_segments$CopyNumber*2

# CNVs bulk
BULK_cnvs <- read.delim(paste(Input_dir,"PMCAHH1-WT/SV/PMCAHH1-WT-C6/purple/circos/PMCAHH1-WT-C6T.cnv.circos", sep =""))
BULK_cnvs <- BULK_cnvs[BULK_cnvs$value < -1 | BULK_cnvs$value > 1,]
BULK_cnvs$CopyNumber <- ifelse(BULK_cnvs$value < -1, "Loss", "Gain")
BULK_cnvs$Sample <-  "CLONE_C6"
colnames(BULK_cnvs)[1] <- "chromosome"
BULK_cnvs[,1] <- gsub(pattern = "hs", replacement = "", x = BULK_cnvs[,1])

# Read the 100kb filtered COBALT files of the PTA samples
readCounts_files <- list(PTA_C6SC1 = paste(Input_dir,"PMCAHH1-WT/PTATO/intermediate/cnvs/readCounts/PMCAHH1-WT-C6SC1.readcounts.filtered.100kb.txt", sep = ""),
                         PTA_C19SC1 = paste(Input_dir,"PMCAHH1-WT/PTATO/intermediate/cnvs/readCounts/PMCAHH1-WT-C19SC1.readcounts.filtered.100kb.txt", sep =""))

readCounts_overview <- data.frame()
for(sample in names(readCounts_files)){
  print(sample)
  print(readCounts_files[[sample]])
  readCounts <- read.delim(readCounts_files[[sample]])
  readCounts$CopyNumber[readCounts$CopyNumber > 4.5] <- 4.5
  readCounts$Sample <- sample
  readCounts_overview <- rbind(readCounts_overview, readCounts)
}
readCounts_overview <- rbind(readCounts_overview, BULK_readcounts100kb)
readCounts_overview$Chromosome <- factor(readCounts_overview$Chromosome, levels = c(1:22, "X", "Y"))

# Read the segment data
readCountsSegments_files <- list(PTA_C6SC1 = paste(Input_dir,"PMCAHH1-WT/PTATO/intermediate/cnvs/readCounts/PMCAHH1-WT-C6SC1.readcounts.segments.txt", sep = ""),
                                 PTA_C19SC1 = paste(Input_dir,"PMCAHH1-WT/PTATO/intermediate/cnvs/readCounts/PMCAHH1-WT-C19SC1.readcounts.segments.txt", sep = ""))

readCountsSegments_overview <- data.frame()
for(sample in names(readCountsSegments_files)){
  print(sample)
  print(readCountsSegments_files[[sample]])
  readCountsSegments <- read.delim(readCountsSegments_files[[sample]])
  #names(readCountsSegments) <- c("Chromosome", "start.pos", "end.pos", "CopyNumber")
  readCountsSegments$mean[readCountsSegments$mean > 4.5] <- 4.5
  readCountsSegments$Sample <- sample
  readCountsSegments_overview <- rbind(readCountsSegments_overview, readCountsSegments)
} 

colnames(readCountsSegments_overview)[colnames(readCountsSegments_overview) == "mean"] <- "CopyNumber"
colnames(readCountsSegments_overview)[colnames(readCountsSegments_overview) == "chrom"] <- "Chromosome"

readCountsSegments_overview_merged <- rbind(readCountsSegments_overview[,c("Chromosome", "start.pos","end.pos","CopyNumber","Sample")], BULK_readcounts_segments[,c("Chromosome", "start.pos","end.pos","CopyNumber","Sample")])

readCountsSegments_overview_merged$Chromosome <- factor(readCountsSegments_overview_merged$Chromosome, levels = c(1:22, "X", "Y"))

# Copy number variants called by PTATO
CNV_files <- list(PTA_C6SC1 = paste(Input_dir, "PMCAHH1-WT/PTATO/intermediate/svs/integration/PMCAHH1/AHH1BULK/PMCAHH1-WT-C6SC1.integrated.cnvs.txt", sep = ""),
                  PTA_C19SC1 = paste(Input_dir, "PMCAHH1-WT/PTATO/intermediate/svs/integration/PMCAHH1/AHH1BULK/PMCAHH1-WT-C19SC1.integrated.cnvs.txt", sep = ""))
# 
CNV_overview <- data.frame()
# 
for(sample in names(CNV_files)){
  print(sample)
  print(CNV_files[[sample]])
  CNVs <- read.delim(CNV_files[[sample]])
  CNVs$Sample <- sample
  CNV_overview <- rbind(CNV_overview, CNVs[,c(1,2,3,6,which(names(CNVs) == "Sample"))])
}


CNV_overview <- rbind(CNV_overview, BULK_cnvs[,c(1:3, 5,6)])
names(CNV_overview) <- c("Chromosome", "start.pos", "end.pos", "CopyNumber", "Sample")
CNV_overview$Chromosome <- factor(CNV_overview$Chromosome, levels = c(1:22, "X", "Y"))

### Color settings
sample_colors <- c( "honeydew4", "darkgrey","honeydew4", "darkgrey")
names(sample_colors) <- names(CNV_files)

CNV_cols <- c("steelblue3","red3","darkorange", "grey60")
names(CNV_cols) <- c("Gain","Loss","cnLOH", "Normal")

### Figure 5a-left panel
p <- ggplot(data = readCounts_overview[readCounts_overview$Chromosome %in% c(1:22),]) +
  geom_rect(data = CNV_overview[CNV_overview$Chromosome %in% c(1:22),], mapping = aes(xmin = start.pos, xmax = end.pos, fill = CopyNumber), ymin = -Inf, ymax = Inf, alpha = 0.25) +
  geom_point(data = readCounts_overview[readCounts_overview$Chromosome %in% c(1:22),], aes(x = (start.pos+end.pos) / 2, y = medianReadCount, col = Sample), alpha = 0.3, size = 0.01, show.legend=FALSE) +
  geom_segment(data = readCountsSegments_overview_merged[readCountsSegments_overview_merged$Chromosome %in% c(1:22),],
               aes(x = start.pos, xend = end.pos, y = CopyNumber, yend = CopyNumber), col = "black") +
  # #geom_hline(yintercept = 0) +
  facet_grid(Sample~Chromosome, space = "free_x", scales = "free", switch="both") +
  theme_bw(base_size = 5) +
  guides(colour = guide_legend(override.aes = list(linewidth=1))) +
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

ggsave(plot = arrangeGrob(q), paste(Output_dir, "FigureS11C.png"), width = 12, height = 7, dpi = 600, units = "cm")
ggsave(plot = arrangeGrob(q), paste(Output_dir, "FigureS11C.pdf"), width = 12, height = 7, units = "cm")


### S11C: Same for IBFM26

# First bin the 1kb COBALT data of the bulk-wgs sample into 100kb windows for fair comparison with the PTA samples
BULK_readcounts1kb <- read.delim(paste(Input_dir, "IBFM26/SV/cobalt/IBFM26-DX2BM-AMLBULKT.cobalt.ratio.tsv", sep = ""))
BULK_readcounts1kb$medianReadcount <- BULK_readcounts1kb$tumorReadCount / mean(BULK_readcounts1kb$tumorReadCount) * 2
BULK_readcounts1kb$CopyNumber <- BULK_readcounts1kb$medianReadcount 

CENTROMERE_FILE	<- paste(Input_dir, "Resources/hg38_centromeres.txt", sep = "")

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

BULK_readcounts1kb_ann <- RemoveCobaltOutliers(BULK_readcounts1kb, CENTROMERES = CENTROMERES_g)
BULK_readcounts1kb_filtered <- BULK_readcounts1kb_ann[BULK_readcounts1kb_ann$FILTER == "PASS",]

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

BULK_readcounts100kb <- MergeBins(BULK_readcounts1kb_filtered , Binsize = 100000, Value_column = 8)
BULK_readcounts100kb$Sample <- "AMLBULK"
BULK_readcounts100kb$CopyNumber[BULK_readcounts100kb$CopyNumber > 4.5] <- 4.5

# Segments bulk

BULK_readcounts_segments <- read.delim(paste(Input_dir,"IBFM26/SV/cobalt/IBFM26-DX2BM-AMLBULKT.cobalt.ratio.pcf", sep = ""))
BULK_readcounts_segments$Sample <-  "AMLBULK"
colnames(BULK_readcounts_segments)[colnames(BULK_readcounts_segments) == "mean"] <- "CopyNumber"
colnames(BULK_readcounts_segments)[colnames(BULK_readcounts_segments) == "chrom"] <- "Chromosome"
BULK_readcounts_segments$CopyNumber <- 2+BULK_readcounts_segments$CopyNumber*2

# CNVs bulk
BULK_cnvs <- read.delim(paste(Input_dir,"IBFM26/SV/purple/circos/IBFM26-DX2BM-AMLBULKT.cnv.circos", sep =""))
BULK_cnvs <- BULK_cnvs[BULK_cnvs$value < -1 | BULK_cnvs$value > 1,]
BULK_cnvs$CopyNumber <- ifelse(BULK_cnvs$value < -1, "Loss", "Gain")
BULK_cnvs$Sample <-  "AMLBULK"
colnames(BULK_cnvs)[1] <- "chromosome"
BULK_cnvs[,1] <- gsub(pattern = "hs", replacement = "", x = BULK_cnvs[,1])

# Read the 100kb filtered COBALT files of the PTA samples
readCounts_files <- list("HSC-PTAP1B8" = paste(Input_dir,"IBFM26/PTATO/intermediate/cnvs/readCounts/IBFM26-DX2BM-HSCPTAP1B8.readcounts.filtered.100kb.txt", sep = ""))

readCounts_overview <- data.frame()
for(sample in names(readCounts_files)){
  print(sample)
  print(readCounts_files[[sample]])
  readCounts <- read.delim(readCounts_files[[sample]])
  readCounts$CopyNumber[readCounts$CopyNumber > 4.5] <- 4.5
  readCounts$Sample <- sample
  readCounts_overview <- rbind(readCounts_overview, readCounts)
}
readCounts_overview$medianReadCount <- readCounts_overview$medianReadCount*2
readCounts_overview <- rbind(readCounts_overview, BULK_readcounts100kb)
readCounts_overview$Chromosome <- factor(readCounts_overview$Chromosome, levels = c(1:22, "X", "Y"))

# Read the segment data
readCountsSegments_files <- list("HSC-PTAP1B8"  = paste(Input_dir,"IBFM26/PTATO/intermediate/cnvs/readCounts/IBFM26-DX2BM-HSCPTAP1B8.readcounts.segments.txt", sep = ""))
readCountsSegments_overview <- data.frame()
for(sample in names(readCountsSegments_files)){
  print(sample)
  print(readCountsSegments_files[[sample]])
  readCountsSegments <- read.delim(readCountsSegments_files[[sample]])
  #names(readCountsSegments) <- c("Chromosome", "start.pos", "end.pos", "CopyNumber")
  readCountsSegments$mean[readCountsSegments$mean > 4.5] <- 4.5
  readCountsSegments$Sample <- sample
  readCountsSegments_overview <- rbind(readCountsSegments_overview, readCountsSegments)
} 

colnames(readCountsSegments_overview)[colnames(readCountsSegments_overview) == "mean"] <- "CopyNumber"
colnames(readCountsSegments_overview)[colnames(readCountsSegments_overview) == "chrom"] <- "Chromosome"

readCountsSegments_overview_merged <- rbind(readCountsSegments_overview[,c("Chromosome", "start.pos","end.pos","CopyNumber","Sample")], BULK_readcounts_segments[,c("Chromosome", "start.pos","end.pos","CopyNumber","Sample")])

readCountsSegments_overview_merged$Chromosome <- factor(readCountsSegments_overview_merged$Chromosome, levels = c(1:22, "X", "Y"))

# Copy number variants called by PTATO
CNV_files <- list("HSC-PTAP1B8" = paste(Input_dir, "IBFM26/PTATO/intermediate/svs/integration/IBFM26/IBFM26-DX2BM-TCELLBULK/IBFM26-DX2BM-HSCPTAP1B8.integrated.cnvs.txt", sep = ""))
# 
CNV_overview <- data.frame()
# 
for(sample in names(CNV_files)){
  print(sample)
  print(CNV_files[[sample]])
  CNVs <- read.delim(CNV_files[[sample]])
  CNVs$Sample <- sample
  CNV_overview <- rbind(CNV_overview, CNVs[,c(1,2,3,6,which(names(CNVs) == "Sample"))])
}

CNV_overview <- rbind(CNV_overview, BULK_cnvs[,c(1:3, 5,6)])
names(CNV_overview) <- c("Chromosome", "start.pos", "end.pos", "CopyNumber", "Sample")
CNV_overview$Chromosome <- factor(CNV_overview$Chromosome, levels = c(1:22, "X", "Y"))

### Color settings
sample_colors <- c( "honeydew4", "darkgrey","honeydew4", "darkgrey")
names(sample_colors) <- names(CNV_files)

CNV_cols <- c("steelblue3","red3","darkorange", "grey60")
names(CNV_cols) <- c("Gain","Loss","cnLOH", "Normal")

### Figure 5a-left panel
p <- ggplot(data = readCounts_overview[readCounts_overview$Chromosome %in% c(1:22),]) +
  geom_rect(data = CNV_overview[CNV_overview$Chromosome %in% c(1:22),], mapping = aes(xmin = start.pos, xmax = end.pos, fill = CopyNumber), ymin = -Inf, ymax = Inf, alpha = 0.25) +
  geom_point(data = readCounts_overview[readCounts_overview$Chromosome %in% c(1:22),], aes(x = (start.pos+end.pos) / 2, y = medianReadCount, col = Sample), alpha = 0.3, size = 0.01, show.legend=FALSE) +
  geom_segment(data = readCountsSegments_overview_merged[readCountsSegments_overview_merged$Chromosome %in% c(1:22),],
               aes(x = start.pos, xend = end.pos, y = CopyNumber, yend = CopyNumber), col = "black") +
  # #geom_hline(yintercept = 0) +
  facet_grid(Sample~Chromosome, space = "free_x", scales = "free", switch="both") +
  theme_bw(base_size = 5) +
  guides(colour = guide_legend(override.aes = list(linewidth=1))) +
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

ggsave(plot = arrangeGrob(q), paste(Output_dir, "FigureS11D.png"), width = 12, height = 6, dpi = 600, units = "cm")
ggsave(plot = arrangeGrob(q), paste(Output_dir, "FigureS11D.pdf"), width = 12, height = 6, units = "cm")

