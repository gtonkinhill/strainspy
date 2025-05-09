#' plot_ani_distribution
#'
#' Generate a Boxplot for ANI Values of a Given Contig and Phenotype
#'
#' This function generates a boxplot to visualize the ANI (Average Nucleotide Identity) values
#' for a specific contig across different levels of a given categorical predictor.
#'
#' @param se SummarizedExperiment. A `SummarizedExperiment` object containing the assay data and metadata.
#' @param genome_file Character. The name of the genome file (contig) to extract from the dataset.
#' @param predictor Character. The phenotype to group samples in the plot.
#' @param plot Logical. If `TRUE` (default), the function returns a ggplot object.
#' If `FALSE`, it returns a data frame containing the plotted data.
#'
#' @return A ggplot object displaying a boxplot of ANI values with a beeswarm overlay across predictor categories.
#'
#' @examples
#' \dontrun{
#' example_meta_path <- system.file("extdata", "example_metadata.csv.gz", package = "strainspy")
#' example_meta <- readr::read_csv(example_meta_path)
#' example_path <- system.file("extdata", "example_sylph_profile.tsv.gz", package = "strainspy")
#' se <- read_sylph(example_path, example_meta)
#' plot_ani_distribution(se, genome_file = 'GCF_025147485.1', predictor = 'Case_status')
#' }
#'
#' @export
plot_ani_distribution <- function(se, genome_file, predictor, plot = TRUE) {

  # Validate input
  if(!inherits(se, 'SummarizedExperiment')) {
    stop('Input must be a SummarizedExperiment object.')
  }

  # Access assay data
  asy <- SummarizedExperiment::assay(se)

  # Access indices of inputted contig and phenotype
  contig <- which(SummarizedExperiment::rowData(se)$Genome_file == genome_file)
  case <- which(colnames(SummarizedExperiment::colData(se)) == predictor)

  # Validate that the contig/phenotype exists
  if(length(contig) == 0 | length(case) == 0) {
    stop('Provide a valid predictor or contig (ex. \'GCF_025147485.1\')')
  }

  # Access the data for the predictor
  case <- SummarizedExperiment::colData(se)[[case]]

  # Create data.frame for plotting
  plot_data <- tibble:tibble(sample = names(asy[contig,]),
                          ani = unname(asy[contig,]),
                          predictor = case)

  # Return data if requested
  if(!plot) {
    return(plot_data)
  }

  custom_colors <- c('#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5')

  custom_colors <- custom_colors[1:length(unique(plot_data$predictor))]

  # Create ggplot with boxplot and beeswarm overlay
  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = predictor, y = ani, fill = predictor)) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.5) +  # semi-transparent boxplot
    ggbeeswarm::geom_quasirandom(width = 0.3, shape = 21, size = 2, stroke = 0.2, alpha = 1) +  # opaque points
    ggplot2::xlab(label = predictor) +
    ggplot2::scale_fill_manual(values = custom_colors) +
    ggplot2::scale_y_continuous(name = 'ANI', limits = c(97, 100)) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = 'none')

  return(plot)

}
