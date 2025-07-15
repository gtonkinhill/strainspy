#' Add or modify metadata in a `SummarizedExperiment` object
#'
#' This function adds or modifies metadata in a `SummarizedExperiment` object.
#' 
#' If `convert = TRUE`, character columns with less than 5 unique (non-NA) values will be automatically
#' converted to factors, with levels ordered alphabetically. For meaningful association testing,
#' however, factor levels should reflect the intended reference level. It is recommended to
#' manually convert such variables to factors with explicitly set levels *before* calling `modify_metadata()`.
#' 
#' @param se A SummarizedExperiment object.
#' @param meta_data A data.frame with sample metadata. **The first column must contain all samples returned by** `colnames(se)`.
#' @param replace Logical. If TRUE, replaces existing colData. Default FALSE.
#' @param convert Logical. If TRUE, attempts to auto-convert types as suited for
#' model fitting (int to numeric, char to factor). Disable if meta_data is fully 
#' configured for model fitting. Default TRUE. See Details.
#'
#'
#' @return Updated SummarizedExperiment object.
#' 
#' @examples
#' \dontrun{
#'   # Add or modify metadata subsequently
#'   example_meta_path <- system.file("extdata", "example_metadata.csv.gz", package = "strainspy")
#'   example_meta <- readr::read_csv(example_meta_path)
#'   example_path <- system.file("extdata", "example_sylph_profile.tsv.gz", package = "strainspy")
#'   se <- read_sylph(example_path)
#'   se <- modify_metadata(se, example_meta)
#' }
#' 
#' @export
modify_metadata <- function(se, meta_data, replace = FALSE, convert = TRUE) {
  # We can't allow duplicate sample names in metadata
  meta_samples <- meta_data[[1]]
  dup_samples <- unique(meta_samples[duplicated(meta_samples)])
  if (length(dup_samples) > 0) {
    n_dups <- length(dup_samples)
    if (n_dups <= 10) {
      stop("Duplicated samples found in meta_data ", paste(dup_samples, collapse = ", "))
    } else {
      stop(
        sprintf("There are %d duplicated samples in meta_data. Showing first 10:\n%s", 
                n_dups, paste(head(dup_samples, 10), collapse = ", "))
      )
    }
  }
  
  # Validate sample names
  se_samples <- colnames(se)
  
  if (length(missing_from_meta <- setdiff(se_samples, meta_samples)) > 0) {
    n_missing <- length(missing_from_meta)
    if (n_missing <= 10) {
      stop("No meta data available for `se` sample(s): ", paste(missing_from_meta, collapse = ", "))
    } else {
      stop(
        sprintf("Meta data is absent for %d 'se' samples. Showing first 10:\n%s", 
                n_missing, paste(head(missing_from_meta, 10), collapse = ", "))
      )
    }
  }
  
  # Covert if requested
  if (convert) {
    meta_data <- auto_coerce_metadata_types(meta_data)
  } else {
    char_cols <- names(meta_data)[sapply(meta_data, is.character)]
    if (length(char_cols) > 0) {
      warning("Character columns detected (", paste(char_cols, collapse = ", "),
              "). Consider converting to factors for modeling.")
    }
  }
  
  # Arrange meta_data in se order and drop sample names
  meta_data_ord <- meta_data[match(se_samples, meta_data[[1]]), , drop = FALSE]
  meta_data_cols <- S4Vectors::DataFrame(meta_data_ord[, -1, drop = FALSE])
  rownames(meta_data_cols) <- se_samples
  
  # Replace colData if requested
  if (replace) {
    SummarizedExperiment::colData(se) <- meta_data_cols
  } else {
    # Extract coldata
    cd = SummarizedExperiment::colData(se)
    # Drop duplicate columns from meta_data_cols
    shared_cols <- intersect(colnames(meta_data_cols), colnames(cd))
    if (length(shared_cols) > 0) {
      warning("The following metadata columns already exist in `se` and will not be modified. To replace with provided meta_data, set run again with `replace = TRUE`:\n\n",
              paste(shared_cols, collapse = ", "))
    }
    
    new_colnames <- setdiff(colnames(meta_data_cols), shared_cols)
    if (length(new_colnames) == 0) {
      warning("No new metadata columns were added, all provided columns already exist in `se`. To replace, run again with `replace = TRUE`.")
    } else {
      new_cols <- meta_data_cols[, setdiff(colnames(meta_data_cols), shared_cols), drop = FALSE]
      merged <- cbind(cd, new_cols)
      
      SummarizedExperiment::colData(se) <- merged
    }
    
  }
  
  return(se)
}

#' Automatically coerce metadata column types for modeling
#'
#' This helper function is intended to standardize metadata columns prior to modeling.
#' It performs the following:
#' - Converts all integer columns to numeric.
#' - Converts character columns to numeric if they contain only numeric-like strings.
#' - Converts character columns to factors if they have fewer than `factor_threshold` unique (non-NA) values.
#' - Leaves all other character columns unchanged.
#' Missing values (NAs) are preserved throughout.
#'
#' @param md A data.frame or tibble of metadata.
#' @param factor_threshold Integer. Maximum number of unique (non-NA) character values to convert to a factor. Default is 5.
#'
#' @return A modified data.frame or tibble with coerced column types.
#'
#' @keywords internal
#' 
#' @importFrom dplyr mutate across where
#' @importFrom S4Vectors na.omit
#' @importFrom magrittr %>%
auto_coerce_metadata_types <- function(md, factor_threshold = 5) {
  md <- md %>%
    dplyr::mutate(
      across(dplyr::where(is.integer), as.numeric),
      across(dplyr::where(is.character), ~ {
        non_na_vals <- .[!is.na(.)]
        suppressWarnings(num_version <- as.numeric(non_na_vals))
        is_numeric_like <- all(!is.na(num_version))
        
        if (is_numeric_like) {
          return(as.numeric(.))  # full column with original NAs preserved
        }
        
        nunique <- length(unique(non_na_vals))
        if (nunique <= factor_threshold) {
          return(factor(., levels = sort(unique(non_na_vals))))
        }
        
        return(.)  # Leave as character
      })
    )
  return(md)
}
