### This script contains functions for plotting mutation patterns

library(tibble)
library(reshape2)
library(grid)

#library(ggh4x) # Could be used a double facet a plot
COLORS6 <- c(
  "#2EBAED", "#000000", "#DE1C14",
  "#D4D2D2", "#ADCC54", "#F0D0CE"
)
COLORS7 <- c(
  "#2EBAED", "#000000", 
  "#E98C7B", "#DE1C14","#D4D2D2", "#ADCC54",
  "#F0D0CE"
)

### This function add text labels  to the stacked bars (eg adds the signature number to the bars), so the signatures are more clear
plot_contribution2 <- function (contribution, signatures = NA, index = NA, coord_flip = FALSE, 
                                mode = c("relative", "absolute"), palette = NA, max_value_to_show = 50) 
{
  mode <- match.arg(mode)
  if (!is.na(index)) {
    contribution <- contribution[, index, drop = FALSE]
  }
  Sample <- Contribution <- Signature <- NULL
  if (mode == "absolute" & !is.na(signatures)) {
    total_signatures <- colSums(signatures)
    abs_contribution <- contribution * total_signatures
  }
  tb <- contribution %>% as.data.frame() %>% tibble::rownames_to_column("Signature") %>% 
    tidyr::pivot_longer(-Signature, names_to = "Sample", 
                        values_to = "Contribution") %>% dplyr::mutate(Sample = factor(Sample, 
                                                                                      levels = unique(Sample)), Signature = factor(Signature, 
                                                                                                                                   levels = unique(Signature)))
  if (mode == "absolute") {
    bar_geom <- geom_bar(stat = "identity", colour = "black")
    text_geom <- geom_text(size = 2, position = position_stack(vjust = 0.5))
    y_lab <- "Absolute contribution \n (no. mutations)"
    y_expansion <- 0.1
    
  }
  else if (mode == "relative") {
    bar_geom <- geom_bar(position = "fill", stat = "identity", 
                         colour = "black")
    text_geom <- geom_text(size = 2,position = position_fill(vjust = 0.5))
    y_lab <- "Relative contribution"
    y_expansion <- 0 
  }
  present_sigs <- tb %>% dplyr::filter(Contribution != 0) %>% 
    dplyr::pull(Signature) %>% unique()
  
  print(tb)
  
  tb$Label <- tb$Signature
  levels(tb$Label) <- c(levels(tb$Label), "")
  tb$Label[tb$Contribution < max_value_to_show] <- ""
  
  # colors are bit different from plot_contribution
  plot <- ggplot(tb[tb$Contribution > 0,], aes(x = Sample, y = Contribution, fill = Signature,
                                               label = gsub(Label, pattern = "SBS", replacement = ""))) + 
    bar_geom  + text_geom +
    labs(x = "", y = y_lab) +
    #scale_y_continuous(expand = expansion(mult = 0), labels = scales::percent) + 
    #scale_fill_discrete(breaks = present_sigs) + 
    scale_y_continuous(expand = expansion(mult = c(0,y_expansion))) +
    theme_bw(base_size = 8) + 
    theme(panel.grid.minor.x = element_blank(), 
          panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), 
          panel.grid.major.y = element_blank(),
          axis.text.x = element_text(angle =45, hjust = 1))
  if (!is.na(palette)) {
    plot <- plot + scale_fill_manual(name = "Signature", 
                                     values = palette)
  }
  if (coord_flip) {
    plot <- plot + coord_flip() + xlim(rev(levels(factor(tb$Sample))))
  }
  else {
    plot <- plot + xlim(levels(factor(tb$Sample)))
  }
  return(plot)
}




## Plot mutation spectra as stacked bar

plot_spectrum_stacked <- function(type_occurrences, mode = "absolute", max_value_to_show = 30, facet = "", total_value = TRUE, round_number = 0, title = ""){
  type_occurrences_stack <- type_occurrences[,c(1,2, 7, 8,4,5,6)]
  ylabel <- "Absolute contribution"
  y_extension <- 0.1
  if(mode == "relative"){
    type_occurrences_stack <- type_occurrences_stack / rowSums(type_occurrences_stack)
    round_number <- 2
    ylabel <- "Relative contribution"
    y_extension <- 0
  }
  
  type_occurrences_stack$Sample <- row.names(type_occurrences_stack)
  type_occurrences_stack$Sample <- factor(type_occurrences_stack$Sample, levels = row.names(type_occurrences_stack))
  if(length(facet > 1)){
    type_occurrences_stack$Facet <- facet
    type_occurrences_stack_m <- melt(type_occurrences_stack, id.vars = c("Sample", "Facet"))
    
  } else {
    type_occurrences_stack_m <- melt(type_occurrences_stack, id.vars = "Sample")
    
  }
  type_occurrences_stack_m$Label <- type_occurrences_stack_m$value
  # The values below max_value_to_show are not plotted
  type_occurrences_stack_m$Label[type_occurrences_stack_m$Label < max_value_to_show] <- ""
  # type_occurrences_total <- aggregate(type_occurrences_stack_m$value ~ type_occurrences_stack_m$Sample, FUN = "sum")
  # colnames(type_occurrences_total) <- c("Sample", "Count")
  if(length(facet) < 2){
    
    plot <- ggplot(data = type_occurrences_stack_m, aes(x = Sample, 
                                                y = value, 
                                                label = as.character(round(as.numeric(Label),round_number )))) + 
      geom_bar(stat = "identity", aes(fill = variable), width = 0.8) +
      geom_text(aes(fill = variable), size = 2, position = position_stack(vjust = 0.5)) +
      #facet_grid(.~Type, space = "free",scales = "free", switch = "both") +
      scale_fill_manual(values = COLORS7) +
      scale_y_continuous(expand = expansion(mult = c(0,y_extension))) +
      theme_classic() +
      ggtitle(title) + 
      labs(x = "Sample", y = ylabel, fill = "Point mutation type")  +
      theme(axis.text.x= element_text(angle = 45, hjust = 1), 
            strip.placement = "outside", 
            strip.background = element_rect(fill = "lightgray"), 
            plot.title = element_text(hjust = 0.5))
    if(total_value == TRUE){
      plot +
        geom_text(aes(label = round(stat(y), round_number), group = Sample), 
                  stat = 'summary', fun = sum, vjust = -1, size = 3) 
    } else {
      plot
    }
  } else {
    ggplot(data = type_occurrences_stack_m, aes(x = Sample, y = value, label = as.character(round(as.numeric(Label),round_number )))) + 
      geom_bar(stat = "identity", aes(fill = variable), width = 0.8) +
      geom_text(aes(fill = variable), size = 2, position = position_stack(vjust = 0.5)) +
      geom_text(aes(label = round(stat(y), 0), group = Sample), 
                stat = 'summary', fun = sum, vjust = -1, size = 3) + 
      facet_grid(.~Facet, scales = "free", space = "free") +
      #facet_grid(.~Type, space = "free",scales = "free", switch = "both") +
      scale_fill_manual(values = COLORS7) +
      scale_y_continuous(expand = expansion(mult = c(0,y_extension))) +
      theme_classic() +
      labs(x = "Sample", y = ylabel, fill = "Point mutation type") +
      ggtitle(title) + 
      theme(axis.text.x= element_text(angle = 45, hjust = 1), 
            strip.placement = "outside", 
            strip.background = element_rect(fill = "lightgray"), 
            plot.title = element_text(hjust = 0.5))
  }

}


plot_96_profile2 <- function (mut_matrix, colors = NA, condensed = FALSE) 
{
  freq <- full_context <- substitution <- context <- NULL
  if (is.na(colors)) {
    colors <- COLORS6
  }
  if (length(colors) != 6) {
    stop("Provide colors vector with length 6", call. = FALSE)
  }
  norm_mut_matrix <- apply(mut_matrix, 2, function(x) x/sum(x))
  tb <- norm_mut_matrix %>% as.data.frame() %>% tibble::rownames_to_column("full_context") %>% 
    dplyr::mutate(substitution = stringr::str_replace(full_context, 
                                                      "\\w\\[(.*)\\]\\w", "\\1"), context = stringr::str_replace(full_context, 
                                                                                                                 "\\[.*\\]", "\\.")) %>% dplyr::select(-full_context) %>% 
    tidyr::pivot_longer(c(-substitution, -context), names_to = "sample", 
                        values_to = "freq") %>% dplyr::mutate(sample = factor(sample, 
                                                                              levels = unique(sample)))
  if (condensed == TRUE) {
    width <- 1
    spacing <- 0
  }
  else {
    width <- 0.6
    spacing <- 0.5
  }
  plot <- ggplot(data = tb, aes(x = context, y = freq, fill = substitution, 
                                width = width)) + geom_bar(stat = "identity", colour = "black", 
                                                           size = 0.2) + scale_fill_manual(values = colors) + facet_grid(sample ~ 
                                                                                                                           substitution) + ylab("Relative contribution") + 
    guides(fill = FALSE) + theme_bw() + theme(axis.title.y = element_text(size = 12, 
                                                                          vjust = 1), axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 12), 
                                              axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5), 
                                              strip.text.x = element_text(size = 9), strip.text.y = element_text(size = 9), 
                                              panel.grid.major.x = element_blank(), panel.spacing.x = unit(spacing, 
                                                                                                           "lines"))
  return(plot)
}


# ggsave doesnt export the facet strip color properly
plot_indel_contexts2 <- function (counts, same_y = FALSE, extra_labels = FALSE, condensed = FALSE, striplabel = "top", stripcolor = "fill",
                                  basesize = 6, legend = FALSE) 
{
  count <- muttype <- muttype_sub <- muttype_total <- sample <- NULL
  counts <- counts %>% as.data.frame() %>% tibble::rownames_to_column("muttype_total") %>% 
    tidyr::separate(muttype_total, c("muttype", "muttype_sub"), 
                    sep = "_(?=[0-9])") %>% dplyr::mutate(muttype = factor(muttype, 
                                                                           levels = unique(muttype))) %>% tidyr::gather(key = "sample", 
                                                                                                                        value = "count", -muttype, -muttype_sub) %>% dplyr::mutate(sample = factor(sample, levels = unique(sample)))
  counts$Label <- ifelse(grepl("deletion", x = counts$muttype) == T, "del", "ins")
  
  # Add the categories to the indel counts (these will used as the facet strip headers)
  indel_categories <- data.frame(muttype = levels(counts$muttype),
                                 category = c(rep("1bp Deletion", 2),
                                              rep("1bp Insertion", 2),
                                              rep(">1bp Deletion at repeats\n(Deletion length)", 4),
                                              rep(">1bp Insertion at repeats\n(Insertion length)", 4),
                                              rep("Deletion with MH\n(Deletion length", 4)),
                                 xaxis = c(rep("Homopolymer length", 2),
                                           rep("Homopolymer length", 2),
                                           rep("Number of repeat units", 4),
                                           rep("Number of repeat units", 4),
                                           rep("Microhomology\nlength", 4)))
  indel_categories$category <- factor(indel_categories$category, levels = c("1bp Deletion", "1bp Insertion",">1bp Deletion at repeats\n(Deletion length)",
                                                                ">1bp Insertion at repeats\n(Insertion length)",
                                                                "Deletion with MH\n(Deletion length"))
  indel_categories$xaxis <- factor(indel_categories$xaxis, levels = c("Homopolymer length", "Number of repeat units",  "Microhomology\nlength"))
  
  counts <- merge(counts,indel_categories, by = "muttype", all.x = T)
  
  nr_muts <- counts %>% dplyr::group_by(sample) %>% dplyr::summarise(nr_muts = round(sum(count)))
  facet_labs_y <- stringr::str_c(nr_muts$sample, "\n(n = ", 
                                 nr_muts$nr_muts, ")")
  names(facet_labs_y) <- nr_muts$sample
  facet_labs_x <- c("C", "T", "C", "T", 2, 3, 4, 
                    "5+", 2, 3, 4, "5+", 2, 3, 4, "5+")
  names(facet_labs_x) <- levels(counts$muttype)
  if (same_y) {
    facet_scale <- "free_x"
  }
  else {
    facet_scale <- "free"
  }
  

  indel_colors <- c( "#FDBE6F", "#FF8001", "#B0DD8B", "#36A12E", "#FDCAB5", 
              "#FC8A6A", "#F14432", "#BC141A", "#D0E1F2", "#94C4DF", 
              "#4A98C9", "#1764AB", "#E2E2EF", "#B6B6D8", "#8683BD", 
              "#61409B")
  strip_colors <-  c(rep("gray90", 8), "#FDBE6F", "#FF8001", "#B0DD8B", "#36A12E", "#FDCAB5", 
                                       "#FC8A6A", "#F14432", "#BC141A", "#D0E1F2", "#94C4DF", 
                                       "#4A98C9", "#1764AB", "#E2E2EF", "#B6B6D8", "#8683BD", 
                                       "#61409B")
  if (extra_labels) {
    title <- stringr::str_c("Deletion           ", "Insertion          ", 
                            "Deletion                                   ", "Insertion                                  ", 
                            "Deletion (MH)")
    x_lab <- stringr::str_c("Homopolymer length                            ", 
                            "Number of repeat units                                                                               ", 
                            "Microhomology length")
  }
  else {
    title <- x_lab <- ""
  }
  if (condensed == TRUE) {
    width <- 1
    spacing <- 0
  }
  else {
    width <- 0.6
    spacing <- 0.5
  }
  
  if(legend == FALSE){
    legend_position <- "none"
  } else {
    legend_position <- "right"
  }
  
  if(striplabel == "top"){
    fig <- ggplot(counts, aes(x = muttype_sub, y = count, fill = muttype, 
                              width = width)) + 
      geom_bar(stat = "identity", col = "white", size = 0.2) + 
      facet_nested(sample ~ xaxis +category + muttype, scales = facet_scale, space = "free_x", labeller = labeller(muttype = facet_labs_x, sample = facet_labs_y) ) + 
      scale_fill_manual(values = indel_colors) + 
      scale_y_continuous(expand = expansion(mult = c(0,0.1))) +   
      theme_classic(base_size = basesize) + 
      labs(fill = "Mutation type", title = title, y = "Indels (#)", x = x_lab) +
      annotation_custom(grid::linesGrob(x = c(0, 0), gp = grid::gpar(lwd = 0.5, col = "lightgrey"))) +
      annotation_custom(grid::linesGrob(y = c(0, 0), gp = grid::gpar(lwd = 1.5, col = "black"))) +
      theme(panel.grid.major.x = element_blank(),
            panel.grid.minor.y = element_blank(), 
            panel.grid = element_blank(),
            panel.spacing.x = unit(spacing, "lines"),
            #panel.border = element_rect(color = "lightgrey", fill = NA),
            strip.background.x = element_rect(color = "white"),
            strip.background.y = element_rect(fill = "gray92", color = "white"),
            axis.text.x = element_text(size = 4),
            strip.placement = "outside", 
            strip.text.y.right = element_text(angle = 0),
            legend.position = legend_position) 

  } else {
    fig <- ggplot(counts, aes(x = muttype_sub, y = count, fill = muttype, 
                              width = width)) + geom_bar(stat = "identity") + 
      facet_nested(sample ~ category +muttype, scales = facet_scale, space = "free_x", labeller = labeller(muttype = facet_labs_x, sample = facet_labs_y), switch = "x") + 
      cale_fill_manual(values = indel_colors) + 
      scale_y_continuous(expand = expansion(mult = c(0,0.1))) +   
      theme_bw(base_size = basesize) + 
      labs(fill = "Mutation type", title = title, y = "Indels (#)", x = x_lab) +
      theme(panel.grid.major.x = element_blank(),
            panel.grid.minor.y = element_blank(), 
            panel.grid = element_blank(),
            panel.spacing.x = unit(spacing, "lines"), 
            strip.background = element_rect(color = "white"),
            strip.placement = "outside", 
            strip.text.y.right = element_text(angle = 0),
            legend.position = legend_position)
  }
  g <- ggplot_gtable(ggplot_build(fig))
  if(striplabel == "top"){
    strip_both <- which(grepl('strip-t', g$layout$name))
  } else {
    strip_both <- which(grepl('strip-b', g$layout$name))
  }
  if(stripcolor == "fill"){
    fills <- strip_colors
  } else{
    fills <- ifelse(grepl("deletion", x = levels(counts$muttype)) == T, "#FB9A99", "#A6CEE3")
  }
  k <- 1
  
  for (i in strip_both) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  plot(g)
  
  #return(fig)
}

# numbers on relative plot don't work for now
plot_main_indel_contexts_stacked <- function (counts, same_y = FALSE, max_value_to_show = 30, round_number = 0, mode = "absolute") 
{
  count <- muttype <- muttype_sub <- muttype_total <- NULL
  counts <- counts %>% as.data.frame() %>% tibble::rownames_to_column("muttype_total") %>% 
    tidyr::separate(muttype_total, c("muttype", "muttype_sub"), 
                    sep = "_(?=[0-9])") %>% dplyr::mutate(muttype = factor(muttype, 
                                                                           levels = unique(muttype)))
  counts_main <- counts %>% dplyr::select(-muttype_sub) %>% 
    dplyr::group_by(muttype) %>% dplyr::summarise_all(list(~sum(.))) %>% 
    tidyr::gather(key = "sample", value = "count", -.data$muttype) %>% 
    dplyr::mutate(sample = factor(sample, levels = unique(sample)))
  nr_muts <- counts_main %>% dplyr::group_by(sample) %>% dplyr::summarise(nr_muts = sum(count))
  facet_labs_y <- stringr::str_c(nr_muts$sample, " (n = ", 
                                 nr_muts$nr_muts, ")")
  names(facet_labs_y) <- nr_muts$sample
  colors <- c("#FDBE6F", "#FF8001", "#B0DD8B", "#36A12E", "#FDCAB5", 
              "#FC8A6A", "#F14432", "#BC141A", "#D0E1F2", "#94C4DF", 
              "#4A98C9", "#1764AB", "#E2E2EF", "#B6B6D8", "#8683BD", 
              "#61409B")
  if (same_y) {
    facet_scale <- "free_x"
  }
  else {
    facet_scale <- "free"
  }
  counts_main$label <- ifelse(counts_main$count > max_value_to_show, counts_main$count, "")
  if(mode == "absolute"){
    fig <- ggplot(counts_main, aes(x = sample, y = count, fill = muttype,  label = as.character(round(as.numeric(label),round_number)))) + 
      geom_bar(stat = "identity",  col = "black", width = 0.75) + 
      geom_text( size = 2, position = position_stack(vjust = 0.5)) +
      geom_text(aes(label = round(stat(y), 0), group = sample), 
                stat = 'summary', fun = sum, vjust = -1, size = 3) +
      labs(x = "", y = "Nr of indels") + 
      scale_fill_manual(values = colors) + 
      scale_y_continuous(expand = expansion(mult = c(0,0.1))) +   
      theme_bw() + 
      theme(panel.grid.major.x = element_blank(), 
            panel.grid.minor.y = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1))
  } else if (mode == "relative") {
    fig <- ggplot(counts_main, aes(x = sample, y = count, fill = muttype,  label = as.character(round(as.numeric(label),round_number)))) + 
      geom_bar(stat = "identity", position = position_fill(), col = "black", width = 0.75) + 
      #geom_text( size = 2, position = position_fill(vjust = 0.5)) +
      #geom_text(aes(label = round(stat(y), 0), group = sample), 
      #          stat = 'summary', fun = sum, vjust = -1, size = 3) +
      labs(x = "", y = "Nr of indels") + 
      scale_fill_manual(values = colors) + 
      scale_y_continuous(expand = expansion(mult = c(0,0.1))) +   
      theme_bw() + 
      theme(panel.grid.major.x = element_blank(), 
            panel.grid.minor.y = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1))
  }

  return(fig)
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


### 96-profiles
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
    "#D4D2D2", "#ADCC54", "#F0D0CE"
  ))
  return(plot_strips)
}

