#' Prepare a caret-compatible prediction data frame from a SummarizedExperiment
#'
#' Formats strain-level data and optional covariates into a data.frame object
#' ready for modeling with \pkg{caret}. The function extracts selected contigs
#' from the `sy` assay matrix, combines them with the specified `outcome` and 
#' optional `covariates` from the metadata, and performs basic validation
#' and cleaning. Metadata is pulled from `SummarizedExperiement::colData(sy)`, 
#' and generally assumes `strainspy::modify_metadata()` was used upstream.
#'
#' @param sy A \link[SummarizedExperiment]{SummarizedExperiment} A 
#' `SummarizedExperiment` object containing the assay data and metadata.
#' @param outcome Character string giving the name of the outcome variable
#'   present in \code{colData(sy)}. This is usually the left-hand side of the
#'   prediction model formula.
#' @param contigs Character vector of contigs (feature) to include as
#'   predictors. Must be present in \code{rownames(sy)}.
#' @param covariates Optional character vector of additional metadata variables
#'   to include as covariates. Default  `NULL` for no covariates.
#' @param use_genome_names Bool. Use `Genome_file` instead of `Contig_name` in 
#' `colnames(output)`. `Genome_file` is extracted from 
#' `SummarizedExperiemnt::elementMetadata(sy)`. If not available, `Contig_name`
#' will be used. Default `True`.
#'
#' @details
#' Character variables are automatically converted to factors.
#' Numeric variables are returned unchanged. For pre processing, refer to 
#' \code{caret::train()}.
#' 
#' @return
#' A \code{data.frame} containing the specified outcome variable, optional
#' covariates, and selected contig features. The object is ready for direct
#' use in \pkg{caret} functions such as \code{caret::train()}. 
#'
#' @examples
#' \dontrun{
#' contigs <- top_hits(fit, coef = 2, method = "BH", alpha = 0.05)$Contig_name
#' df <- prep_for_prediction(
#'   sy = sy,
#'   outcome = "Case_status",
#'   contigs = contigs,
#'   covariates = c("Sex", "BMI")
#' )
#' }
#'
#' @export
prep_for_prediction <- function(sy, outcome, contigs, covariates = NULL, use_genome_names = T) {
  meta <- as.data.frame(SummarizedExperiment::colData(sy))
  sy_mx <- SummarizedExperiment::assay(sy)
  
  if (!outcome %in% names(meta)) stop("Outcome variable '", outcome, "' not found in metadata.")
  
  if (!all(contigs %in% rownames(sy_mx))) stop("Some contigs not found in assay data.")
  
  sy_mx <- sy_mx[rownames(sy_mx) %in% contigs, , drop = FALSE]
  sy_mx <- t(as.matrix(sy_mx))
  
  if(use_genome_names){
    contig_genomes = as.data.frame(SummarizedExperiment::elementMetadata(sy))
    genomes = contig_genomes$Genome_file[match(colnames(sy_mx), contig_genomes$Contig_name)]
    if(any(duplicated(genomes))) stop("Found duplicated genomes.")
    colnames(sy_mx) = make.names(genomes)
  }
  
  df <- data.frame(outcome = meta[[outcome]])
  names(df)[1] <- outcome
  
  if (!is.null(covariates)) {
    missing_cov <- setdiff(covariates, names(meta))
    if (length(missing_cov) > 0)
      stop("Missing covariates in metadata: ", paste(missing_cov, collapse = ", "))
    df <- cbind(df, meta[, covariates, drop = FALSE])
  }
  
  df <- cbind(df, as.data.frame(sy_mx))

  char_cols <- sapply(df, is.character)
  if (any(char_cols)) {
    df[char_cols] <- lapply(df[char_cols], as.factor)
    message("Converted ", sum(char_cols), " character variable(s) to factors.")
  }
  
  if (any(!is.finite(as.matrix(df[, sapply(df, is.numeric), drop = FALSE])))) {
    warning("Detected Inf or NaN values in some of the numeric columns. Check before feeding to <CARET> or other prediction pipeline.")
  }
  
  if (any(is.na(df))) {
    warning("Detected missing (NA) values in the prepared data frame.")
  }

  if (is.character(df[[outcome]]) || is.factor(df[[outcome]])) df[[outcome]] <- as.factor(df[[outcome]])
  
  message("Prepared data: ", nrow(df), " samples and ", ncol(df) - 1, " predictors.")
  return(df)
}
