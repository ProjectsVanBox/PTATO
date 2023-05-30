### Supplemental Figure S12: copy number and DAF plots FA patients

#library(VariantAnnotation)
library(StructuralVariantAnnotation)
#library(stringr)
library(ggplot2)

input_dir <- "/path/to/MendeleyData/directory/"
output_dir <- "/path/to/output/directory/"

if(dir.exists(output_dir) == F){dir.create(output_dir)}

Metadata <- read.delim(paste(input_dir, "Table_S1.txt", sep = ""))
Metadata_FA <- Metadata[Metadata$Phenotype == "FanconiAnemia" & Metadata$Germline == FALSE,]

#' Simple SV type classifier
simpleEventType <- function(gr) {
  pgr = partner(gr)
  return(ifelse(seqnames(gr) != seqnames(pgr), "CTX", # inter-chromosomosal
                ifelse(strand(gr) == strand(pgr), "INV",
                       ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS", 
                              ifelse(xor(start(gr) < start(pgr), strand(gr) == "-"), "DEL",
                                     "DUP")))))
}


# Coverage needs to contain chromosome - start - end - value
# DAF needs to contain chromosome - start - end - value

Combined_ReadCount_DAF_plot <- function(Coverage, 
                                        DAF, 
                                        Chromosomes = c(1:22,"X", "Y"),
                                        PlotColors = c("red", "blue"),
                                        Max_copynumber = 4, Title = ""){
  
  Coverage[,1] <- factor(Coverage[,1], levels = Chromosomes)
  colnames(Coverage)[4] <- "Value"
  Coverage$Type <- "CopyNumber"
  Coverage$Value[Coverage$Value > Max_copynumber] <- Max_copynumber
  
  DAF[,1] <- factor(DAF[,1], levels = Chromosomes)
  colnames(DAF)[4] <- "Value"
  DAF$Type <- "DAF"
  
  Merged_data <- rbind(Coverage, DAF)
  colnames(Merged_data) <- c("Chromosome", "Start", "End", "Value", "Type")
  Merged_data$Color <- "Even"
  Merged_data$Color[Merged_data$Chromosome %in% c(seq(1,21, by = 2), "Y")] <- "Uneven"
  
  # Dummy is necessary to keep all y-axiss the same for each different facet (from 0-4 for copynumber and 0-0.5 for DAF)
  Dummy <- data.frame(Chromosome = rep(1, 4),
                      Start = rep(1,4), 
                      End = rep(2,4),
                      Value = c(0, Max_copynumber,0, 0.5), 
                      Type = c("CopyNumber","CopyNumber", "DAF","DAF"),
                      Color = "Uneven")
  
  Plot <- ggplot(Merged_data, aes(x = (Start+End) / 2, y = Value, color = Color)) +
    geom_point(size = 0.02, alpha = 0.8) +
    geom_blank(data=Dummy, aes(x = (Start+End) / 2, y = Value, color = Color)) +
    facet_grid(Type~Chromosome, scales = "free", space = "free_x", switch = "both") + 
    scale_color_manual(values = c("black", "dimgray")) +
    #coord_cartesian(ylim = c(0,4)) +
    theme_minimal(base_size = 6) +
    theme(axis.text.x = element_blank(),
          panel.grid.major.y = element_line(linewidth  = 0.1),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background.y = element_rect(fill = "lightgray", colour  = "white"), 
          strip.placement = "outside",
          strip.text.x = element_text(size = 4),
          panel.spacing.x  = unit(0, "lines"), 
          panel.background = element_rect(fill = "gray97", colour = NA),
          panel.border = element_rect(fill = NA, colour = "white"),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none") +
    labs(x = "Genomic position (bp)", y = "")
  
  if(Title != ""){
    Plot <- Plot + ggtitle(Sample) +
      theme(plot.title = element_text(hjust = 0.5, size = 5))
  }
  
  return(Plot)
}

### Read the coverage and DAF data
for(Sample in Metadata_FA$Sample){
  print(Sample)
  
  Coverage_file <- list.files(path = paste(input_dir, Metadata_FA$Individual[Metadata_FA$Sample == Sample], sep = ""), 
                              pattern = paste(Sample,".readcounts.filtered.1mb.txt", sep = ""), recursive = T, full.names = T)
  print(Coverage_file)
  
  DAF_file <- list.files(path = paste(input_dir, Metadata_FA$Individual[Metadata_FA$Sample == Sample], sep = ""), 
                              pattern = paste(Sample,".baf.binned.1mb.txt", sep = ""), recursive = T, full.names = T)
  print(DAF_file)
  
  if(length(Coverage_file) >0  & length(DAF_file> 0)){
    Coverage_sample <- read.delim(Coverage_file)
    
    DAF_sample <- read.delim(DAF_file)
    
    Combined_ReadCount_DAF_plot(Coverage = Coverage_sample[,c(1:3,5)], DAF = DAF_sample[,c(1:4)], Max_copynumber = 5, Title = Sample)
    ggsave(paste(output_dir, "CN_", Sample, ".pdf",sep = ""), width = 8.5, height = 3.5, units = "cm")
    ggsave(paste(output_dir, "CN_", Sample, ".png",sep = ""), width = 8.5, height = 3.5, units = "cm", dpi = 300)
    
  } else {
    print("!!! FILES NOT FOUND")
  }
}

Combined_ReadCount_DAF_plot(Coverage = Coverage_sample[,c(1:3,5)], DAF = DAF_sample[,c(1:4)], Max_copynumber = 5)

