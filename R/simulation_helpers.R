#' Add or modify metadata in an `SummarizedExperiment` object
#'
#' This function adds or modifies metadata in `SummarizedExperiment` object read using functions such as `read_sylph()` and `read_metaphlan()`.
#'
#' @param se SummarizedExperiment object generated from `strainspy` read functions.
#' @param meta_data data.frame. A tibble or data frame containing sample metadata.
#' @param replace bool. If T, the meta data in se (i.e., `se@coldata`) will be replaced with the provided meta_data file. Default F.
#' The **first column must contain sample names** that match exactly with the `colnames(se)`.
#' @return SummmarizedExperiment object `se`, updated with the new meta data.
modify_metadata <- function(se, meta_data, replace = F) {
  meta_samples <- unique(meta_data[[1]])
  if (length(missing_from_meta <- setdiff(se@colData$Sample_file, meta_samples)) > 0) {
    stop("The following samples from 'se' are not in 'meta_data': ", paste(missing_from_meta, collapse = ", "))
  }
  
  if (length(missing_from_se <- setdiff(meta_samples, se@colData$Sample_file)) > 0) {
    stop("The following samples from 'meta_data' are not in 'se': ", paste(missing_from_se, collapse = ", "))
  }
  
  # check and reorder as necessary
  asy = SummarizedExperiment::assay(se)
  if ( !all(names(asy[1,]) == se@colData@rownames) ) {
    stop("Mismatch between se rownames and assay order")
  }
  
  if(replace == T){
    se@colData@listData = list(se@colData@listData$Sample_file)
  }
  
  # we can't have the same colname in se@colData and meta_data
  test_cols = which(colnames(meta_data) %in% colnames(se@colData))
  if(length(test_cols) >0) {
    warning(paste("Columns:", colnames(meta_data)[test_cols], "exist(s) in se@colData and will be dropped without merging. Set replace = T to replace instead" ))
    meta_data = meta_data[,-test_cols]
  }
  
  if(ncol(meta_data) > 0) {
    new_metadata <- S4Vectors::DataFrame(base::merge(se@colData, meta_data,
                                                     by.x = "Sample_file",
                                                     by.y = names(meta_data)[1],
                                                     all.x = TRUE,
                                                     sort = FALSE)) # Keeps all Sample_file entries from se
    
    se@colData@listData = as.list(new_metadata)
  } else {
    warning("No new metadata to add, returning the same se object.")
  }
  
  return(se)
}

#' Update the model fit by calling the specified fit function only a subset of features. Same design will be fitted.
#'
#' @param fit Output fit from `glmZiBFit()`, `glmFit()` or `caseControlFit()`
#' @param update_idx Vector of indices that need to be updated
#' @param se SummarizedExperiment. A `SummarizedExperiment` object containing the assay data and metadata.
#' @param scale_continous Binary specifying whether to rescale numeric values. Default T
#' @param min_identity Only for `caseControlFit()`. A numeric value specifying the minimum identity threshold to consider (default=0.98).
#' @param method Character. The method to use for fitting the model. Either 'glmmTMB' (default) or 'gamlss'. Only applicable for `glmZiBFit()`.
#' @param family A `glmmTMB` family object. Defaults to `glmmTMB::ordbeta()`. Only applicable for `glmFit()`
#' @param nthreads An integer specifying the number of (CPUs or workers) to use. Defaults to 1.
#' @param BPPARAM Optional `BiocParallelParam` object. If not provided, the function will configure an appropriate backend automatically.
#'        
#' @return Updated version of the fit object 
#'
#' @import SummarizedExperiment
#' @importFrom glmmTMB glmmTMB
update_fit <- function(fit, update_idx, 
                       se, scale_continous=TRUE,
                       min_identity=0.98, method='glmmTMB', family=glmmTMB::ordbeta(), 
                       nthreads=1,  BPPARAM=NULL) {
  
  if (!inherits(se, "SummarizedExperiment")) {
    stop("`se` must be a SummarizedExperiment object.")
  }
  
  if (!inherits(fit, "betaGLM")) {
    stop("`fit` must be a betaGLM object.")
  }
  
  # what was fitted?
  fit_type = "Unknown"
  fits = c("glmZiBFit", "glmFit", "caseControlFit")
  mtch = which(fits %in% as.character(fit@call[[1]]))
  if(length(mtch) != 0) fit_type = fits[mtch]
  
  if(fit_type == "Unknown") {
    stop("`fit` type must be one of glmZiBFit, glmFit or caseControlFit")
  }
  
  
  # modify data to subset update_idx
  se_subset = se[update_idx,]
  design = as.formula(fit@call[["design"]])
  
  switch(fit_type,
         "glmZiBFit" = {
           cat("Updating", length(update_idx), "values by calling glmZiBFit\n")
           fit_u <- glmZiBFit(se_subset, design = design, nthreads = nthreads, scale_continous = scale_continous, BPPARAM = BPPARAM)
         },
         "glmFit" = {
           cat("Updating", length(update_idx), "values by calling glmFit\n")
           fit_u <- glmFit(se_subset, design = design, nthreads = nthreads, scale_continous = scale_continous, BPPARAM = BPPARAM, family = family)
         },
         "caseControlFit" = {
           cat("Updating", length(update_idx), "values by calling caseControlFit\n")
           fit_u <- caseControlFit(se = se_subset, min_identity = min_identity, design = design, nthreads = nthreads, scale_continous = scale_continous, BPPARAM = BPPARAM)
         },
         {
           # Default case (optional)
           stop("Invalid fit_type")
         }
  )
  
  # let's update fit
  if ("coefficients" %in% slotNames(fit)) {
    fit@coefficients[update_idx, ] = fit_u@coefficients
  }
  
  if ("std_errors" %in% slotNames(fit)) {
    fit@std_errors[update_idx, ] = fit_u@std_errors
  }
  
  if ("p_values" %in% slotNames(fit)) {
    fit@p_values[update_idx, ] = fit_u@p_values
  }
  
  if ("residuals" %in% slotNames(fit)) {
    fit@residuals[update_idx, ] = fit_u@residuals
  }
  
  # These are not present in all models, check NULLs as well
  if ("zi_coefficients" %in% slotNames(fit)) {
    if(!is.null(fit@zi_coefficients)){
      fit@zi_coefficients[update_idx, ] = fit_u@zi_coefficients
    }
  }
  
  if ("zi_std_errors" %in% slotNames(fit)) {
    if(!is.null(fit@zi_std_errors)){
      fit@zi_std_errors[update_idx, ] = fit_u@zi_std_errors
    }
  }
  
  if ("zi_p_values" %in% slotNames(fit)) {
    if(!is.null(fit@zi_p_values)){
      fit@zi_p_values[update_idx, ] = fit_u@zi_p_values
    }
  }
  
  
  
  return(fit)
  
}


#' Generate a confusion matrix for a fitted model
#'
#' @param top_hits Output from `tophits(fit)`
#' @param gt_contigs Vector of ground truth contigs 
#' @param all_contigs Vector of all contigs
#' @param print_cm logical. Print the confusion matrix (default = TRUE)
#'        
#' @return Confusion Matrix.
get_confusion_mx <- function(top_hits, gt_contigs, all_contigs, print_cm = TRUE) {
  cm = matrix(c( sum(top_hits$Contig_name %in% gt_contigs),
                 sum( !(top_hits$Contig_name %in% gt_contigs) ), 
                 sum( !(gt_contigs %in% top_hits$Contig_name) ), 
                 sum( !(all_contigs %in% c(top_hits$Contig_name, gt_contigs) ) ) ), 
              nrow = 2, 
              byrow = TRUE,
              dimnames = list("Detected" = c("Positive", "Negative"), 
                              "Actual" = c("Positive", "Negative")))
  if(print_cm == TRUE) {
    print(cm)
  }
  
  return(cm)
}

#' Generate a Manhattan plot from a betaGLM object, with the option to add ground truth 
#'
#' This function is copied from plot_manhattan() 
#'
#' @param object A `betaGLM` object.
#' @param coef The number of the coefficient from which to generate the plot (default=2).
#' @param plot If set to false, the function will return a tibble with data used to generate the plot.
#' @param ground_truth A vector of ground truth contigs
#' @return A `ggplot` object showing the Manhattan plot.
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate
plot_manhattan_gt <- function(object, coef=2, taxonomy=NULL, method = "HMP",
                              alpha=0.05, levels = c("Phylum", "Genus", "Species", "Strain"),
                              colour_level = "Phylum", plot=TRUE, ground_truth = NULL) {
  
  # Validate input
  if (!inherits(object, "betaGLM")) {
    stop("Input must be a betaGLM object.")
  }
  
  # Create tidy tibbles for plotting
  if (!is.null(taxonomy)) {
    # Add taxonomic informed adjusted p-values
    plot_data <- hadjust(object, coef=coef, taxonomy = taxonomy,
                         method = method, index_range = TRUE)[, c('Level', 'Model', 'Name', 'p_adjust', 'index_min', 'index_max')]
    
    tax_levels <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
    plot_data$Level[!plot_data$Level %in% tax_levels] <- "Strain"
    plot_data$Level <- factor(plot_data$Level, levels=tax_levels)
  } else {
    # Add sequence level adjusted p-values
    hits <- top_hits(object, coef=coef, method = method, alpha = 1)
    
    plot_data <- hits[,c(colnames(object@row_data), "p_adjust")] |>
      tibble::add_column(Model="Beta", .before=1)
    
    if ('zi_p_adjust' %in% colnames(hits)) {
      plot_data <- rbind(plot_data,
                         hits[,c(colnames(object@row_data), "zi_p_adjust")] |>
                           dplyr::rename(p_adjust = "zi_p_adjust") |>
                           tibble::add_column(Model="Zero-Inflated", .before=1))
    }
    
    colnames(plot_data)[[2]] <- "Name"
    plot_data <- plot_data |> tibble::add_column(Level="Sequence", .before=1)
    
    plot_data$index <- 1:nrow(plot_data)
  }
  
  plot_data$log_p_adjust <- -log10(plot_data$p_adjust)
  
  if(!is.null(ground_truth)) plot_data$ground_truth = as.factor(plot_data$Name %in% ground_truth)
  
  # Create the Manhattan plot
  if (!is.null(taxonomy)){
    # add colours for chosen level
    ranges <- plot_data |>
      dplyr::filter(Level==colour_level) |>
      dplyr::group_by(Name) |>
      dplyr::summarise(
        min_index = min(index_min),
        max_index = max(index_max)
      ) |>
      dplyr::rename(Phylum = Name)
    
    plot_data <- plot_data |>
      dplyr::cross_join(ranges) |>  # Cross join
      dplyr::filter(index_min >= min_index & index_min <= max_index) |>
      dplyr::select(-min_index, -max_index)  # Remove unnecessary columns
    
    plot_data <- plot_data |> dplyr::filter(Level %in% levels)
    
    # Return data if requested
    if (!plot){
      return(plot_data)
    }
    
    # Define your custom color palette
    custom_colors <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#b15928',
                       '#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5')
    
    custom_colors <- custom_colors[1:length(unique(plot_data$Phylum))]
    
    gg <- ggplot2::ggplot(plot_data, aes(x = Name, y = log_p_adjust, colour=Phylum)) +
      ggplot2::geom_segment(aes(y = log_p_adjust, x = index_min - 0.5, xend = index_max + 0.5), size=5) +
      ggplot2::facet_grid(Level ~ Model, scales = "free_x") +
      geom_hline(aes(yintercept = -log10(alpha)), col='red', linetype = "dotted") +
      scale_color_manual(values = custom_colors) +
      ggthemes::theme_clean(base_size = 16) +
      ggplot2::theme(plot.background = element_blank()) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(face = "bold"),
      ) +
      ggplot2::ylab("-log10(p-value)") +
      ggplot2::xlab("Strain")
    
  } else {
    
    # Return data if requested
    if (!plot){
      return(plot_data)
    }
    
    ## Only edited this function for now
    
    if(!is.null(ground_truth)){
      
      gg <- ggplot2::ggplot(plot_data, aes(x = Name, y = log_p_adjust, col = ground_truth)) +
        ggplot2::geom_point() +
        ggplot2::facet_wrap( ~ Model) +
        geom_hline(aes(yintercept = -log10(alpha)), col='red', linetype = "dotted") +
        ggthemes::theme_clean(base_size = 16) +
        ggplot2::theme(plot.background = element_blank()) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          strip.text = ggplot2::element_text(face = "bold"),
        ) +
        ggplot2::ylab("-log10(p-value)") +
        ggplot2::xlab("Strain")
      
    } else {
      
      gg <- ggplot2::ggplot(plot_data, aes(x = Name, y = log_p_adjust)) +
        ggplot2::geom_point() +
        ggplot2::facet_wrap( ~ Model) +
        geom_hline(aes(yintercept = -log10(alpha)), col='red', linetype = "dotted") +
        ggthemes::theme_clean(base_size = 16) +
        ggplot2::theme(plot.background = element_blank()) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          strip.text = ggplot2::element_text(face = "bold"),
        ) +
        ggplot2::ylab("-log10(p-value)") +
        ggplot2::xlab("Strain")
      
    }
    
  }
  
  return(gg)
}

#' Generate a Violin plot for a subset of contigs for a given categorical phenotype 
#'
#' @param se  SummarizedExperiment. A `SummarizedExperiment` object containing the assay data and metadata.
#' @param phenotype One variable in `colnames(se@colData)`, except `Sample_File`. The selected variable must be categorical
#' @param contigs A vector of contigs to generate violin plots
#' @param contig_names A vector of names for the contigs to use in the plot, useful when the original contig names are too long. Must be the same length and order as `contig_names`. Default NULL.
#' @param drop_ANI_zeros Bool. If TRUE, ANI 0 values will not be included in the violin plot. Default FALSE. A warning will be generated to inform that ANI 0 rows are dropped.
#' @param plot If set to false, the function will return a tibble with data used to generate the plot.
plot_violin <- function(se, phenotype, contigs, contig_names = NULL, drop_ANI_zeros = FALSE, plot = TRUE){
  if(length(phenotype) != 1) {
    stop("Only one phenotype can be plotted at a time")
  }
  
  if(! (phenotype %in% colnames(se@colData)) ){
    stop("Phenotype not found in colnames(se@colData)")
  }

  if( ! is.factor(se@colData[[phenotype]])) {
    stop("Phenotype must be a factor")
  }
  
  chk_contigs = contigs %in% rownames(se)
  if( ! all(chk_contigs)  ){
    cat(paste(contigs[!chk_contigs], '\n', sep = ""))
    stop("Above contigs are missing from rownames(se)")
  }
  
  tmp = as.data.frame(t(as.matrix(SummarizedExperiment::assay(se[contigs, ]))))
  tmp = cbind(se@colData[[phenotype]], tmp)
  
  if(!is.null(contig_names)){
    colnames(tmp) = c(phenotype, contig_names)
  } else {
    colnames(tmp) = c(phenotype, colnames(tmp)[-1])
  }
  
  df_long <- tidyr::pivot_longer(tmp, cols = -phenotype, names_to = "Contig", values_to = "ANI")
  
  if(drop_ANI_zeros){
    df_long$ANI[df_long$ANI == 0] = NA
  }
  
  phenotype <- rlang::sym(phenotype)
  
  if(!plot){
    return(df_long)
  } else {
    ggplot(df_long, aes(x = Contig, y = ANI, fill = !!phenotype)) +
      geom_violin(trim = FALSE, alpha = 0.6) +
      theme_minimal() +
      labs(x = "Contigs", y = "ANI") +
      scale_fill_brewer(palette = "Set2") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  
}

#TODO: Add Nikhil's box plot code here









