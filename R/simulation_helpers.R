#' Add or modify metadata in an `SummarizedExperiment` object
#'
#' This function adds or modifies metadata in `SummarizedExperiment` object read using functions such as `read_sylph()` and `read_metaphlan()`.
#'
#' @param se SummarizedExperiment object generated from `strainspy` read functions.
#' @param meta_data data.frame. A tibble or data frame containing sample metadata.
#' The **first column must contain sample names** that match exactly with the `colnames(se)`.
#' @return SummmarizedExperiment object `se`, updated with the new meta data.
modify_metadata <- function(se, meta_data) {
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
  
  new_metadata <- S4Vectors::DataFrame(base::merge(se@colData, meta_data,
                                                   by.x = "Sample_file",
                                                   by.y = names(meta_data)[1],
                                                   all.x = TRUE,
                                                   sort = FALSE)) # Keeps all Sample_file entries from se
  
  se@colData@listData = as.list(new_metadata)
  
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
