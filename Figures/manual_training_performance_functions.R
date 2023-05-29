# Function to read the training data
read_train_data = function(fname, label){
    df <- readRDS(fname)
    if (nrow(df)){
        df$val <- label
    }
    return(df)
}

# Function to subset the classes of a df to be of equal size
subset_classes_equal = function(df){
    n_per_class = min(table(df$val))
    df_subset = df %>%
        dplyr::group_by(val) %>% 
        dplyr::slice_sample(n = n_per_class) %>% 
        dplyr::ungroup()
    return(df_subset)
}

# Plot the importance of a random forest
plot_importance = function(rf){
    imp = importance(rf) %>% 
        as.data.frame() %>% 
        rownames_to_column()
    imp_fig = ggplot(imp, aes(x = reorder(rowname, MeanDecreaseGini), y = MeanDecreaseGini)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        labs(x = "Feature") +
        theme_classic()
    return(imp_fig)
}

# Function to determine the performance of the rf on subsets of the data
subset_rf = function(df, omit){
    n = nrow(df)
    n_v = seq(100, n, by = 20)
    n_v = rep(n_v, 3)
    cost_tb = purrr::map(n_v, subset_rf_i, df, omit, n) %>%
        bind_rows() %>%
        dplyr::group_by(n) %>%
        dplyr::summarise(across(everything(), mean))
    return(cost_tb)
}

# Helper function of subset_rf
subset_rf_i = function(n_sub, df, omit, n){
    
    # Subset df
    df_sub = df[sample(n, n_sub),]
    
    # Fit the model
    if (omit | !sum(is.na(df_sub))){
        rf = randomForest(as.factor(val) ~ ., data=df_sub, mtry=4)
    } else{
        train_data_imputed = rfImpute(val ~ ., data = df_sub)
        rf = randomForest(as.factor(val) ~ ., data=train_data_imputed, mtry=4)
    }
    
    # Determine cost
    conf = rf$confusion
    TP = conf["noPTA", "noPTA"]
    FP = conf["PTA", "noPTA"]
    TN = conf["PTA", "PTA"]
    FN = conf["noPTA", "PTA"]
    oob = rf$err.rate[rf$ntree,1][[1]]
    cost = data.frame(n = n_sub, oob = oob, TP = TP, FP = FP, TN = TN, FN = FN) %>%
        dplyr::mutate(precision = TP / (TP + FP),
                      precision = ifelse(is.nan(precision), 1, precision),
                      recall = TP / (TP + FN),
                      specificity = TN / (TN + FP),
                      accuracy = (TP + TN) / (TP + FP + FN + TN),
                      F1 = 2 * precision * recall / (precision + recall),)
    return(cost)
}

# Determine the accuracy of the rf when using different cutoffs
determine_cutoff_accuracy = function(p, res){
    TP <- nrow(res[res$val == "noPTA" & res$pta_prob <= p,])
    FP <- nrow(res[res$val == "PTA" & res$pta_prob <= p,])
    FN <- nrow(res[res$val == "noPTA" & res$pta_prob > p,])
    TN <- nrow(res[res$val == "PTA" & res$pta_prob > p,])
    precision <- TP / (TP + FP)
    recall <- TP / (TP + FN)
    specificity <- TN / (TN + FP)
    accuracy <- (TP + TN) / (TP + FP + FN + TN)
    F1 = 2 * (precision * recall) / (precision + recall)
    noPTA_class_error = FN / (TP + FN)
    PTA_class_error = FP / (TN + FP)
    res = data.frame("cutoff" = p, TP, FP, FN, TN, precision, recall, specificity, accuracy, F1, noPTA_class_error, PTA_class_error)
    return(res)
}

#Function to calculate the area under the curve of the rf roc curve.
calculate_auc = function(tb){
    tb = dplyr::select(tb, precision, recall)
    max_x = max(tb$precision)
    curve = approxfun(tb, rule = 2, ties = max)
    auc = integrate(curve, 0, max_x)$value
    return(auc)
}

# Perform mutational pattern analyses on a mutation matrix
do_pta_mutpatterns = function(mut_mat, signatures, group){
    
    
    # Look at profile
    prof_fig = plot_96_profile(mut_mat, condensed = TRUE)
    
    # Look at signature contribution
    fit_res <- fit_to_signatures(mut_mat, signatures)
    contri_fig <- plot_contribution(fit_res$contribution, coord_flip = FALSE, mode = "absolute")
    
    fit_res2 <- fit_to_signatures(mut_mat, signatures[,c("PTA", "HSPC")])
    contri_fig2 <- plot_contribution(fit_res2$contribution, coord_flip = FALSE, mode = "absolute")
    
    strict_refit <- fit_to_signatures_strict(mut_mat, signatures, max_delta = 0.05)
    fit_res_strict <- strict_refit$fit_res
    contri_strict_fig <- plot_contribution(fit_res_strict$contribution,
                                           coord_flip = FALSE,
                                           mode = "absolute")
    
    # Look at cosine similarity
    cos_sim_samples_signatures <- cos_sim_matrix(mut_mat, signatures)
    cossim_fig <- plot_cosine_heatmap(cos_sim_samples_signatures,
                                      cluster_rows = FALSE, cluster_cols = FALSE)
    
    saveRDS(list(mut_mat, fit_res, fit_res2, strict_refit, cos_sim_samples_signatures), paste0(group, "_mutpatterns.rds"))
    pdf(paste0(group, "_mutpatterns.pdf"))
    print(prof_fig)
    print(contri_fig)
    print(contri_fig2)
    print(contri_strict_fig)
    print(cossim_fig)
    dev.off()
}

# Perform mutational pattern analyses on a mutation matrix
do_pta_mutpatterns_indel = function(indel_counts, signatures, group){
    
    # Look at profile
    prof_fig = plot_indel_contexts(indel_counts, condensed = TRUE)
    
    # Look at signature contribution
    fit_res <- fit_to_signatures(indel_counts, signatures)
    contri_fig <- plot_contribution(fit_res$contribution, coord_flip = FALSE, mode = "absolute")
    
    # fit_res2 <- fit_to_signatures(indel_counts, signatures[,c("PTA", "HSPC")])
    # contri_fig2 <- plot_contribution(fit_res2$contribution, coord_flip = FALSE, mode = "absolute")
    
    strict_refit <- fit_to_signatures_strict(indel_counts, signatures, max_delta = 0.05)
    fit_res_strict <- strict_refit$fit_res
    contri_strict_fig <- plot_contribution(fit_res_strict$contribution,
                                           coord_flip = FALSE,
                                           mode = "absolute")
    
    # Look at cosine similarity
    cos_sim_samples_signatures <- cos_sim_matrix(indel_counts, signatures)
    cossim_fig <- plot_cosine_heatmap(cos_sim_samples_signatures,
                                      cluster_rows = FALSE, cluster_cols = FALSE)
    
    saveRDS(list(indel_counts, fit_res, strict_refit, cos_sim_samples_signatures), paste0(group, "_indel_mutpatterns.rds"))
    pdf(paste0(group, "_mutpatterns.pdf"))
    print(prof_fig)
    print(contri_fig)
    print(contri_strict_fig)
    print(cossim_fig)
    dev.off()
}

# Function to get complement of base
my_complement = function(base){
    base_comp = case_when(base == "T" ~ "A",
                  base == "A" ~ "T",
                  base == "C" ~ "G",
                  base == "G" ~ "C")
        #str_replace("T", "A") %>% 
        #str_replace("C", "G")
    return(base_comp)
}

# Helper function of count_indel_contexts. Modified here to work without the ref/alt columns
.count_indel_contexts_gr = function(gr, categories) {
    
    # Check gr is not empty
    if (length(gr) == 0) {
        categories <- categories %>%
            dplyr::mutate(count = 0) %>%
            dplyr::select(-muttype, -muttype_sub)
        return(categories)
    }
    
    # Check context has previously been set.
    gr_colnames <- colnames(mcols(gr))
    if (!all(c("muttype", "muttype_sub") %in% gr_colnames)) {
        stop("The GRanges object does not contain the columns `muttype`` and `muttype_sub`.
             Did you forget to run `get_indel_context`?", call. = FALSE)
    }
    
    
    # Classify the number of repeat units/ homopolymer length / microhomology length
    # to either 5+ or 6+ depending on whether the indel is a ins or del.
    id_context <- dplyr::tibble("muttype" = gr$muttype, "muttype_sub" = gr$muttype_sub) %>%
        dplyr::mutate(
            muttype_sub = ifelse(muttype_sub >= 6, "6+", muttype_sub),
            muttype_sub = ifelse(grepl("insertion|microhomology", muttype) & muttype_sub >= 5,
                                 "5+", muttype_sub
            ),
            muttype_sub = as.character(muttype_sub)
        ) # Ensures column type for later joining
    
    
    id_context_count <- id_context %>%
        dplyr::group_by(muttype, muttype_sub) %>%
        dplyr::summarise(count = dplyr::n())
    id_context_count_full <- dplyr::left_join(categories,
                                              id_context_count,
                                              by = c("muttype", "muttype_sub")
    ) %>%
        dplyr::select(-muttype, -muttype_sub)
    # colnames(id_context_count_full)[3] = name
    return(id_context_count_full)
}

# Only the helper function is modified. This function is only supplied here, 
# so it calls the helper function from the global env, instead of the regular one.
count_indel_contexts_mod = function(vcf_list) {
    
    # These variables use non standard evaluation.
    # To avoid R CMD check complaints we initialize them to NULL.
    muttype <- muttype_sub <- NULL
    
    categories <- INDEL_CATEGORIES
    
    # Turn grl into list if needed.
    if (inherits(vcf_list, "CompressedGRangesList")) {
        vcf_list <- as.list(vcf_list)
    }
    
    # Count contexts per sample
    if (inherits(vcf_list, "list")) {
        counts_l <- purrr::map(vcf_list, .count_indel_contexts_gr, categories)
        counts <- do.call(cbind, counts_l)
        colnames(counts) <- names(vcf_list)
    } else if (inherits(vcf_list, "GRanges")) {
        counts <- .count_indel_contexts_gr(vcf_list, categories)
        colnames(counts) <- "My_sample"
    } else {
        .not_gr_or_grl(vcf_list)
    }
    counts <- cbind(categories, counts)
    counts[is.na(counts)] <- 0
    counts <- counts %>%
        tidyr::unite("muttype_total", muttype, muttype_sub) %>%
        tibble::column_to_rownames("muttype_total") %>%
        as.matrix()
    return(counts)
}

plot_conf = function(conf){
    conf_fig = ggplot(conf, aes(x = predicted, y = val, fill = n)) +
        geom_raster() +
        geom_text(aes(label = n), size = 2) +
        scale_fill_distiller(palette = "Blues", limits = c(0, 1), direction = 1) +
        labs(x = "Predicted", y = "Real value", fill = "Ratio") +
        my_theme
    return(conf_fig)
}

plot_conf_donor = function(donor_conf){
    donor_conf_fig = ggplot(donor_conf, aes(x = predicted, y = val, fill = n)) +
        geom_raster() +
        geom_text(aes(label = n)) +
        facet_grid(DONOR_ID ~ .) +
        theme_classic()
    return(donor_conf_fig)
}

my_theme = theme(text = element_text(size = 6, family = "Arial"),
                 legend.background = element_rect(fill="transparent", colour=NA),
                 panel.background = element_rect(fill='transparent'),
                 legend.key = element_rect(fill="transparent", colour=NA),
                 axis.text = element_text(size = rel(1), 
                                          colour = "black"),
                 axis.ticks.x = element_line(colour = "black", size = 0.25),
                 axis.ticks.y = element_line(colour = "black", size = 0.25),
                 axis.line = element_line(size = 0.25),
                 strip.background = element_blank(),
                 plot.title = element_text(hjust = 0.5),
                 axis.ticks.length=unit(0.05, "cm"),
                 legend.key.size = unit(0.3, 'cm'),
                 plot.margin = unit(c(5.5,5.5,5.5,5.5), "pt"),
                 legend.text=element_text(size=6),
                 legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
)

spectrum_theme = theme(text = element_text(size = 6, family = "Arial"),
                       legend.background = element_rect(fill="transparent", colour=NA),
                       legend.key = element_rect(fill="transparent", colour=NA),
                       axis.text = element_text(size = rel(1), 
                                                colour = "black"),
                       axis.ticks.x = element_line(colour = "black", size = 0.25),
                       axis.ticks.y = element_line(colour = "black", size = 0.25),
                       axis.ticks.length=unit(0.05, "cm"),
                       axis.line = element_line(size = 0.25),
                       plot.title = element_text(hjust = 0.5),
                       strip.background = element_rect(size = 0.5),
                       panel.border = element_rect(size = 0.25),
                       legend.key.size = unit(0.3, 'cm'),
                       plot.margin = unit(c(5.5,5.5,5.5,5.5), "pt"),
                       legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
                       axis.title.y = element_text(size = 6, vjust = 1),
                       axis.text.y = element_text(size = 6),
                       axis.title.x = element_text(size = 6),
                       axis.text.x = element_text(size = 2.5, angle = 90, vjust = 0.5),
                       strip.text.x = element_text(size = 6),
                       strip.text.y = element_text(size = 6),
                       legend.text=element_text(size=6),
                       panel.spacing.y = unit(0.2, "lines")
)

plot_indel_contexts <- function(counts, same_y = FALSE, extra_labels = FALSE, condensed = FALSE, label_nr_muts = TRUE) {
    # These variables use non standard evaluation.
    # To avoid R CMD check complaints we initialize them to NULL.
    count <- muttype <- muttype_sub <- muttype_total <- sample <- NULL
    
    # Separate muttype and muttype_sub. Then make data long
    counts <- counts %>%
        as.data.frame() %>%
        tibble::rownames_to_column("muttype_total") %>%
        tidyr::separate(muttype_total, c("muttype", "muttype_sub"), sep = "_(?=[0-9])") %>%
        dplyr::mutate(muttype = factor(muttype, levels = unique(muttype))) %>%
        tidyr::gather(key = "sample", value = "count", -muttype, -muttype_sub) %>% 
        dplyr::mutate(sample = factor(sample, levels = unique(sample)))
    
    # Count nr mutations. (This is used for the facets)
    nr_muts <- counts %>%
        dplyr::group_by(sample) %>%
        dplyr::summarise(nr_muts = round(sum(count)))
    
    # Create facet texts
    if (label_nr_muts){
        facet_labs_y <- stringr::str_c(nr_muts$sample, " (n = ", nr_muts$nr_muts, ")")
    } else{
        facet_labs_y <- nr_muts$sample
    }
    names(facet_labs_y) <- nr_muts$sample
    facet_labs_x <- c("1: C", "1: T", "1: C", "1: T", 2, 3, 4, "5+", 2, 3, 4, "5+", 2, 3, 4, "5+")
    names(facet_labs_x) <- levels(counts$muttype)
    
    # Set plotting parameters
    if (same_y) {
        facet_scale <- "free_x"
    } else {
        facet_scale <- "free"
    }
    
    # Add optional extra labels
    if (extra_labels) {
        title <- stringr::str_c(
            "Deletion           ",
            "Insertion          ",
            "Deletion                                   ",
            "Insertion                                  ",
            "Deletion (MH)"
        )
        x_lab <- stringr::str_c(
            "Homopolymer length                            ",
            "Number of repeat units                                                                               ",
            "Microhomology length"
        )
    } else {
        title <- x_lab <- ""
    }
    
    # Change plotting parameters based on whether plot should be condensed.
    if (condensed == TRUE) {
        width <- 1
        spacing <- 0
    } else {
        width <- 0.6
        spacing <- 0.5
    }
    
    # Create figure
    fig <- ggplot(counts, aes(x = muttype_sub, y = count, fill = muttype, width = width)) +
        geom_bar(stat = "identity") +
        facet_grid(sample ~ muttype,
                   scales = facet_scale, space = "free_x",
                   labeller = labeller(muttype = facet_labs_x, sample = facet_labs_y)
        ) +
        scale_fill_manual(values = INDEL_COLORS) +
        theme_bw() +
        labs(fill = "Mutation type", title = title, y = "Nr of indels", x = x_lab) +
        theme(
            panel.grid.major.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.spacing.x = unit(spacing, "lines")
        )
    
    return(fig)
}
