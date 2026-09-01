#' Strainspy Wrapper
#'
#' A wrapper to run a full strainspy workflow with explicit input validation.
#'
#' @param meta_path Path to the metadata file (CSV format).
#' @param sylph_path Path to the sylph output file (TSV format).
#' @param design_formula One-sided formula (or string coercible to formula) for
#' the design matrix, e.g. `~ Condition + Age + (1|Subject)`.
#' @param coef Coefficient index to report in `top_hits()`. This determines
#' which design term is tested.
#' @param min_nonzero Minimum number of samples for strain presence filtering.
#' If `NULL` (default), this is set to `ceiling(0.10 * n_samples)`.
#' @param taxonomy_path Optional path to a taxonomy file (TSV format).
#' If provided, taxonomy columns are appended to top hits.
#' @param output_path Optional path to save top hits as CSV. If `NULL`, no file is written.
#' @param alpha Significance threshold used by `top_hits()`.
#' @param p_adjust_method Multiple-testing correction method for `top_hits()`.
#' @param model_method Backend used by `glmZiBFit()`: `"glmmTMB"` or `"gamlss"`.
#' @param scale_continuous If `TRUE`, continuous metadata covariates are z-scored
#' by `glmZiBFit()`.
#' @param nthreads Number of threads for parallel processing (default: all available cores).
#' @param return_vcov If `TRUE`, include per-feature variance-covariance matrices in fit.
#' @param attach_taxonomy If `TRUE` and `taxonomy_path` is provided, append taxonomy columns to top hits.
#' @param overwrite If `FALSE`, abort when `output_path` already exists.
#' @return A data frame of top hits. If `taxonomy_path` is supplied and
#' `attach_taxonomy = TRUE`, taxonomy columns are appended.
#' @examples
#' \donttest{
#' meta_path <- system.file("extdata", "example_metadata.csv.gz", 
#' package = "strainspy")
#' sylph_path <- system.file("extdata", "example_sylph_profile.tsv.gz", 
#' package = "strainspy")
#' taxonomy_path <- system.file("extdata", "example_taxonomy.tsv.gz", 
#' package = "strainspy")
#' 
#' # Decision 1: choose predictors for the design formula
#' # Decision 2: choose which coefficient index to report in top_hits()
#' # For ~ Case_status + Age_at_collection, coef = 2 corresponds to Case_status
#' 
#' out <- strainspy(
#'   meta_path = meta_path,
#'   sylph_path = sylph_path,
#'   design_formula = "~ Case_status + Age_at_collection",
#'   coef = 2,
#'   taxonomy_path = taxonomy_path,
#'   output_path = tempfile(fileext = ".csv"),
#'   nthreads = 2,
#'   alpha = 1
#' )
#' head(out)
#' }
#' @export
strainspy <- function(meta_path, sylph_path, design_formula, coef,
                      min_nonzero = NULL,
                      taxonomy_path = NULL, output_path = NULL,
                      alpha = 0.05,
                      p_adjust_method = "holm",
                      model_method = "glmmTMB",
                      scale_continuous = TRUE,
                      nthreads = 2,
                      return_vcov = FALSE,
                      attach_taxonomy = TRUE,
                      overwrite = TRUE) {
  if (!is.character(meta_path) || length(meta_path) != 1 || !nzchar(meta_path)) {
    stop("`meta_path` must be a non-empty character path to a metadata CSV file.")
  }
  if (!file.exists(meta_path)) {
    stop("Metadata file not found at: ", meta_path)
  }

  if (!is.character(sylph_path) || length(sylph_path) != 1 || !nzchar(sylph_path)) {
    stop("`sylph_path` must be a non-empty character path to a sylph TSV file.")
  }
  if (!file.exists(sylph_path)) {
    stop("Sylph profile file not found at: ", sylph_path)
  }

  if (missing(design_formula)) {
    stop("`design_formula` is required (e.g., `~ Condition + Age`).")
  }

  if (missing(coef)) {
    stop("`coef` is required. Provide the coefficient index to extract from `top_hits()`.")
  }

  if (!is.null(taxonomy_path)) {
    if (!is.character(taxonomy_path) || length(taxonomy_path) != 1 || !nzchar(taxonomy_path)) {
      stop("`taxonomy_path` must be NULL or a non-empty character path.")
    }
    if (!file.exists(taxonomy_path)) {
      stop("Taxonomy file not found at: ", taxonomy_path)
    }
  }

  if (!is.null(min_nonzero)) {
    if (!is.numeric(min_nonzero) || length(min_nonzero) != 1 || is.na(min_nonzero) || min_nonzero < 1) {
      stop("`min_nonzero` must be NULL or a single positive integer.")
    }
    min_nonzero <- as.integer(min_nonzero)
  }

  if (!is.numeric(coef) || length(coef) != 1 || is.na(coef) || coef < 1 || coef %% 1 != 0) {
    stop("`coef` must be a single positive integer index.")
  }
  coef <- as.integer(coef)

  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) || alpha < 0 || alpha > 1) {
    stop("`alpha` must be a single numeric value between 0 and 1.")
  }

  if (!is.numeric(nthreads) || length(nthreads) != 1 || is.na(nthreads) || nthreads < 1) {
    stop("`nthreads` must be a single positive integer.")
  }
  nthreads <- as.integer(nthreads)

  if (!is.logical(attach_taxonomy) || length(attach_taxonomy) != 1 || is.na(attach_taxonomy)) {
    stop("`attach_taxonomy` must be TRUE or FALSE.")
  }
  if (!is.logical(scale_continuous) || length(scale_continuous) != 1 || is.na(scale_continuous)) {
    stop("`scale_continuous` must be TRUE or FALSE.")
  }
  if (!is.logical(return_vcov) || length(return_vcov) != 1 || is.na(return_vcov)) {
    stop("`return_vcov` must be TRUE or FALSE.")
  }
  if (!is.logical(overwrite) || length(overwrite) != 1 || is.na(overwrite)) {
    stop("`overwrite` must be TRUE or FALSE.")
  }

  design <- tryCatch({
    if (inherits(design_formula, "formula")) {
      design_formula
    } else if (is.character(design_formula) && length(design_formula) == 1) {
      stats::as.formula(design_formula)
    } else {
      stop("`design_formula` must be a one-sided formula or a single formula string.")
    }
  }, error = function(e) {
    stop("Failed to parse `design_formula`: ", conditionMessage(e))
  })

  if (length(design) != 2) {
    stop("`design_formula` must be one-sided (e.g., `~ Condition + Age`). Do not include a response variable.")
  }

  design_vars <- all.vars(design)
  if (length(design_vars) == 0) {
    stop("`design_formula` does not contain any predictor terms.")
  }

  meta_data <- tryCatch({
    readr::read_csv(meta_path, show_col_types = FALSE)
  }, error = function(e) {
    stop("Failed to read metadata file `", meta_path, "`: ", conditionMessage(e))
  })

  if (!is.data.frame(meta_data) || ncol(meta_data) < 2) {
    stop("`meta_data` must contain at least two columns: sample ID in column 1 and at least one metadata variable.")
  }

  missing_design_cols <- setdiff(design_vars, colnames(meta_data))
  if (length(missing_design_cols) > 0) {
    stop(
      "Metadata is missing predictor column(s) required by `design_formula`: ",
      paste(missing_design_cols, collapse = ", "),
      ".\nAvailable metadata columns are: ",
      paste(colnames(meta_data), collapse = ", "),
      "."
    )
  }

  se <- read_sylph(sylph_path)

  se <- tryCatch({
    modify_metadata(se, meta_data, replace = TRUE, convert = TRUE)
  }, error = function(e) {
    stop(
      "Failed to attach metadata to `SummarizedExperiment`: ", conditionMessage(e),
      "\nThe first column of metadata must contain all sample names from the sylph file."
    )
  })

  if (is.null(min_nonzero)) {
    min_nonzero <- as.integer(max(1L, ceiling(0.10 * ncol(se))))
  }

  missing_in_coldata <- setdiff(design_vars, colnames(SummarizedExperiment::colData(se)))
  if (length(missing_in_coldata) > 0) {
    stop(
      "After metadata merge, design variable(s) are missing from `colData(se)`: ",
      paste(missing_in_coldata, collapse = ", "),
      "."
    )
  }

  se <- filter_by_presence(se, min_nonzero = min_nonzero)
  if (nrow(se) == 0) {
    stop("No strains remain after `filter_by_presence()` with min_nonzero = ", min_nonzero, ".")
  }

  fit <- glmZiBFit(
    se,
    design,
    nthreads = nthreads,
    scale_continuous = scale_continuous,
    method = model_method,
    return_vcov = return_vcov
  )

  top_hits_df <- top_hits(fit, coef = coef, method = p_adjust_method, alpha = alpha)

  if (!is.null(taxonomy_path) && isTRUE(attach_taxonomy) && nrow(top_hits_df) > 0) {
    taxonomy <- read_taxonomy(taxonomy_path)
    taxonomy_df <- as.data.frame(taxonomy)
    idx <- match(top_hits_df$Genome_file, taxonomy$Genome)
    tax_cols <- setdiff(colnames(taxonomy_df), c("Genome", "index"))
    taxonomy_subset <- taxonomy_df[idx, tax_cols, drop = FALSE]
    top_hits_df <- cbind(top_hits_df, taxonomy_subset)
  }

  if (!is.null(output_path)) {
    if (!is.character(output_path) || length(output_path) != 1 || !nzchar(output_path)) {
      stop("`output_path` must be NULL or a non-empty character path.")
    }
    out_dir <- dirname(output_path)
    if (!dir.exists(out_dir)) {
      stop("Output directory does not exist: ", out_dir)
    }
    if (file.exists(output_path) && !isTRUE(overwrite)) {
      stop("Output file already exists and `overwrite = FALSE`: ", output_path)
    }
    readr::write_csv(top_hits_df, output_path)
  }

  top_hits_df
}