#' ani_distance
#'
#' This function calculates the distance between samples in a `SummarizedExperiment` object representing strain level ANI data.
#'
#' @param se A `SummarizedExperiment` object.
#' @param taxonomy An optional `data.frame` object containing the taxonomy information for the samples in `se`.
#' The result of running read_taxonomy on the taxonomy file associated with the ANI data.
#' If provided, the ANI values will be collapsed at the specified taxonomic level before calculating the distance matrix. (default=NULL)
#' @param tax_level A character string specifying the taxonomic level at which to collapse the ANI values. (default="Species")
#' @param collapse_method A character string specifying the method to use to collapse the ANI values. Must be one of 'representative', 'mean' or 'max'. (default='representative')
#'
#' @return A distance matrix.
#'
#' @examples
#' \dontrun{
#' example_path <- system.file("extdata", "example_sylph_profile.tsv.gz", package = "strainspy")
#' tax_path <- system.file("extdata", "example_taxonomy.tsv.gz", package = "strainspy")
#' se <- read_sylph(example_path)
#' tax <- read_taxonomy(tax_path)
#' distance_matrix <- ani_distance(se)
#' }
#'
#' @export
ani_distance <- function(se, taxonomy=NULL, tax_level="Species", collapse_method='representative') {

  # Check that the input is a SummarizedExperiment object
  if (!inherits(se, "SummarizedExperiment")) {
    stop("`se` must be a SummarizedExperiment object.")
  }

  # Collapse at specified taxa level is taxonomy is provided
  if (!is.null(taxonomy)){
    fullm <- as.matrix(SummarizedExperiment::assays(se)[[1]])
    grouping <- tax[[tax_level]][match(rowData(se)$Genome_file, tax$Genome)]
    if (collapse_method=='representative') {
    ani_vals <- t(do.call(rbind, purrr::map(unique(grouping), ~ {
      index <- which.max(rowMeans(fullm[grouping == .x, , drop=FALSE]))
      return(fullm[index, , drop=FALSE])
    })))
    } else if (collapse_method=='mean') {
      ani_vals <- t(do.call(rbind, purrr::map(unique(grouping), ~ {
        return(colMeans(fullm[grouping == .x, , drop=FALSE]))
      })))
    } else if (collapse_method=='max') {
      ani_vals <- t(do.call(rbind, purrr::map(unique(grouping), ~ {
        return(colMaxs(fullm[grouping == .x, , drop=FALSE]))
      })))
    } else {
      stop("Invalid collapse method. Must be 'representative', 'mean' or 'max'.")
    }
  } else {
    ani_vals <- t(as.matrix(SummarizedExperiment::assays(se)[[1]]))
  }

  # Calculate distances
  d <- dist(ani_vals, method = "euclidean")

  # Return the distance matrix
  return(d)
}
