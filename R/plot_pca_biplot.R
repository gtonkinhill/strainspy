#' plot_pca_biplot
#'
#' Perform PCA and Generate a Biplot
#'
#' This function performs Principal Component Analysis (PCA) on a SummarizedExperiment object
#' and generates a PCA biplot colored by a specified phenotype.
#'
#' @param se A `SummarizedExperiment` object containing expression data.
#' @param coef A character string or integer specifying the column in the meta data by which the plot is coloured
#' @param plot A logical value indicating whether to return the PCA plot (`TRUE`) or PCA data (`FALSE`). Default is `TRUE`.
#'
#' @return If `plot = TRUE`, returns a `ggplot2` PCA plot. If `plot = FALSE`, returns a tibble with PCA results.
#'
#' @examples
#' \dontrun{
#' example_meta_path <- system.file("extdata", "example_metadata.csv.gz", package = "strainspy")
#' example_meta <- readr::read_csv(example_meta_path)
#' example_path <- system.file("extdata", "example_sylph_profile.tsv.gz", package = "strainspy")
#' se <- read_sylph(example_path, example_meta)
#' plot_pca_biplot(se, coef = 3)
#' }
#'
#' @export
plot_pca_biplot <- function(se, coef = 2, plot = TRUE) {

  asy<-SummarizedExperiment::assay(se)

  # PCA
  PC <- prcomp(t(as.matrix(asy)), scale = TRUE)

  # Calculate the percent variance from the first two PC
  pve <- PC$sdev^2/sum(PC$sdev^2)
  xLab <- paste0("PC1 (",round(pve[1]*100,2),"%)")
  yLab <- paste0("PC2 (",round(pve[2]*100,2),"%)")

  #Format PC data
  if(is.null(rownames(PC$x))){
    rownames(PC$x)<-1:nrow(PC$x)
  }
  dat <- tibble::as_tibble(PC$x)

  # Access the phenotype by which the plot is coloured
  col <- SummarizedExperiment::colData(se)[, coef]
  dat$pred <- col

  # Return plot data if requested
  if(!plot) {
    return(dat)
  }

  # Create a PCA plot
  plot <- ggplot2::ggplot(
    data = dat,
    mapping = ggplot2::aes(x=PC1, y=PC2, color = pred)
    ) +
    ggplot2::geom_point()+
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme_classic()+
    ggplot2::labs(x=xLab,y=yLab)

  # Add a gradient colour if phenotype is continuous
  if(is.numeric(col)) {
    plot <- plot + ggplot2::scale_color_viridis_c(option = 'plasma')
  }

  return(plot)
}
