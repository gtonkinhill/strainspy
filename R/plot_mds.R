#' Plot MDS of ANI data
#'
#' Performs classical multidimensional scaling (MDS) on a `SummarizedExperiment` object
#' and plots the result, optionally colored by a phenotype.
#'
#' @param se A `SummarizedExperiment` object containing ANI data.
#' @param coef Integer. Index of the column in `colData(se)` to use as the grouping or color variable (default is 2).
#' @param taxonomy Optional taxonomy table to use in distance calculation.
#' @param palette A vector of color hex codes for numeric phenotype gradients. Defaults to viridis palette.
#' @param plot Logical. If `TRUE`, returns a `ggplot` object. If `FALSE`, returns the underlying MDS coordinates as a tibble.
#'
#' @return A `ggplot` object showing the MDS plot, or a tibble of MDS coordinates if `plot = FALSE`.
#'
#' @examples
#' \dontrun{
#' example_meta_path <- system.file("extdata", "example_metadata.csv.gz", package = "strainspy")
#' example_meta <- readr::read_csv(example_meta_path)
#' example_path <- system.file("extdata", "example_sylph_profile.tsv.gz", package = "strainspy")
#' se <- read_sylph(example_path, example_meta)
#' plot_mds(se, coef = 2)
#' }
#'
#' @export
plot_mds <- function(se,
                     coef = 2,
                     taxonomy = NULL,
                     palette = c("#DAA520", "#3CB371", "#483D8B"),
                     plot = TRUE) {

  # Access phenotype metadata
  col <- SummarizedExperiment::colData(se)[,coef]

  # Grab the distance matrix and perform MDS
  mds <- ani_distance(se = se, taxonomy = taxonomy) |>
    stats::cmdscale() |>
    tibble::as_tibble() |>
    # Add the phenotype data for colouring
    dplyr::mutate(phenotype = col, Dim1 = V1, Dim2 = V2) |>
    dplyr::select(-V1,-V2)

  if(!plot) {
    return(mds)
  }

  # MDS Plot
  plot <- ggplot2::ggplot(mds, mapping = ggplot2::aes(
    x = Dim1, y = Dim2, colour = phenotype)) +
    ggplot2::geom_point() +
    ggthemes::theme_few()

  # Gradient colouring if continuous data
  if(is.numeric(col)) {
    plot <- plot + scale_colour_gradientn(colours = palette)
  }

  return(plot)
}
