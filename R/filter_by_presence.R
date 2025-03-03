#' filter_by_presence
#'
#' This function filters a SummarizedExperiment object to keep only rows (features/strains)
#' that are present in a specified number of samples.
#' The same rows are filtered from all assays in the SummarizedExperiment.
#'
#' @param se A `SummarizedExperiment` object.
#' @param min_nonzero Integer. The minimum number of non-zero entries required to retain a row. Default is 10.
#'
#' @return A filtered `SummarizedExperiment` object with only the rows that meet the criteria.
#'
#' @import SummarizedExperiment
#'
#' @examples
#' \dontrun{
#' example_path <- system.file("extdata", "example_sylph_profile.tsv.gz", package = "strainspy")
#' se <- read_sylph(example_path)
#'
#' # Filter to keep only rows with at least 10 non-zero entries in the Adjusted_ANI assay
#' filtered_se <- filter_by_presence(se)
#'
#' # View the filtered object
#' filtered_se
#' }
#'
#' @export
filter_by_presence <- function(se, min_nonzero = 10) {

  # Check that the input is a SummarizedExperiment object
  if (!inherits(se, "SummarizedExperiment")) {
    stop("`se` must be a SummarizedExperiment object.")
  }

  # Isn't it convenient to use an abundance threshold instead (like: 0.05 * ncol)?
  if(min_nonzero %% 1 != 0) {
    stop("`min_zero` must be an integer.")
  }

  # Count the number of non-zero entries in each row of the assay
  nonzero_counts <- Matrix::rowSums(SummarizedExperiment::assays(se)[[1]] != 0)

  # Identify which rows have at least 'min_nonzero' non-zero entries
  rows_to_keep <- nonzero_counts >= min_nonzero

  # Filter the assays, rowData, and colData in the SummarizedExperiment
  cat("Retained", sum(rows_to_keep), "rows after filtering\n")
  filtered_se <- se[rows_to_keep, ]

  # Return the filtered SummarizedExperiment
  return(filtered_se)
}
