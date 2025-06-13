#' Generate a Volcano Plot from a betaGLM Object
#'
#' This function creates a volcano plot from the results of a betaGLM analysis,
#' highlighting the relationship between effect size (coefficients) and statistical significance (-log10 p-values).
#'
#' @importFrom viridis viridis
#' 
#' @param fit A betaGLM object containing model results.
#' @param coef Integer, specifying the phenotype index of the betaGLM to extract p values and effect size for plotting. Default is 2.
#' @param alpha Float indicating signficance level
#' @param palette A vector of user defined palette to colour plot. Default is the viridis palette.
#' @param label Logical, whether to label most signifcant points of the plot with Species name.
#' @param contig_names Optional character vector to replace long contig names in the plot.
#' Must match the length and order of `top_hits(fit, coef = coef, alpha = alpha)$Contig_name`.
#' If NULL (default), shortening will be attempted automatically.
#' @param plot Logical, whether to return a ggplot object (default: TRUE). If FALSE, returns the processed data.
#'
#'
#' @return A ggplot2 object if plot = TRUE, otherwise a dataframe containing the relevant statistics.
#'
#' @examples
#' \dontrun{
#' example_meta_path <- system.file("extdata", "example_metadata.csv.gz", package = "strainspy")
#' example_meta <- readr::read_csv(example_meta_path)
#' example_path <- system.file("extdata", "example_sylph_profile.tsv.gz", package = "strainspy")
#' se <- read_sylph(example_path, example_meta)
#' se <- filter_by_presence(se)
#' design <- as.formula(" ~ Case_status + Sex")
#' fit <- glmZiBFit(se, design, nthreads=parallel::detectCores())
#' plot_volcano(fit_glmmtmb, coef = 2, label = T)
#' }
#'
#' @export
plot_volcano <- function(fit, coef = 2, alpha = 0.5,
                         palette = rev(viridis::viridis(3)),
                         label = FALSE, contig_names = NULL,
                         plot = TRUE) {
  
  # Validate input
  if(!inherits(fit, 'betaGLM')) {
    stop('Input must be a betaGLM object.')
  }
  
  # Access p values and effect sizes
  hits <- top_hits(fit, coef = coef, alpha = 1)
  
  # contig names provided
  if(!is.null(contig_names)){
    n_contigs = nrow(top_hits(fit, coef = coef, alpha = alpha))  
    if(length(contig_names) != nrow(top_hits(fit, coef = coef, alpha = alpha))){
      stop(sprintf("Provided %d `contig_names`, but there are %d top_hits contigs at alpha = %.3f", 
                   length(contig_names), n_contigs, alpha))
    }
  }
  
  # Checks if multiple models exist in the fit
  if ('zi_p_adjust' %in% colnames(hits)) {
    
    # Creating a tibble that can be facetted by Model
    p <- hits |>
      tidyr::pivot_longer(cols = c(p_value, zi_p_value), names_to = 'Model', values_to = 'p') |>
      dplyr::select(Genome_file, Model, p) |>
      dplyr::mutate(Model = ifelse(stringr::str_detect(Model, 'zi'), 'ZiB', 'Beta'))
    p_adj <- hits |>
      tidyr::pivot_longer(cols = c(p_adjust, zi_p_adjust), names_to = 'Model', values_to = 'p_adj') |>
      dplyr::select(Genome_file, Model, p_adj) |>
      dplyr::mutate(Model = ifelse(stringr::str_detect(Model, 'zi'), 'ZiB', 'Beta'))
    hits <- hits |>
      tidyr::pivot_longer(cols = c(coefficient,zi_coefficient),
                          names_to = 'Model', values_to = 'Coefficient') |>
      dplyr::mutate(Model = ifelse(stringr::str_detect(Model, 'zi'), 'ZiB', 'Beta')) |>
      dplyr::inner_join(p, by = c('Genome_file', 'Model')) |>
      dplyr::inner_join(p_adj, by = c('Genome_file', 'Model'))
    
    # Return the data if requested
    if(!plot) {
      return(hits)
      
    }
    
    # Plot facetted volcano plots
    plot <- hits |> ggplot2::ggplot(ggplot2::aes(x = Coefficient, y = -log10(p), colour = p_adj)) +
      ggplot2::geom_point() +
      ggplot2::facet_grid(~ Model, scale = 'free_x') +
      ggthemes::theme_few() +
      ggplot2::scale_color_gradientn(colors = palette) +
      ggplot2::theme(panel.spacing = ggplot2::unit(1, 'cm')) +
      ggplot2::xlab('Effect Size') +
      ggplot2::ylab('-log10(p)')
    
    if(label) {
      plot <- plot + ggrepel::geom_text_repel(data = hits[hits$p_adj < alpha,],
                                              mapping = ggplot2::aes(label = clean_contig_names(Contig_name)),
                                              colour = 'black',
                                              size = 2.5,
                                              alpha = 0.8,
                                              max.overlaps = Inf)
    }
  } else { # Plot one-model volcano plot
    hits <- hits |>
      dplyr::mutate(log10p = log10(p_value)) |>
      dplyr::mutate(sig = ifelse(p_adjust < alpha, T, F))
    
    if(!plot) {
      return(hits)
    }
    
    plot <- hits |> ggplot2::ggplot(ggplot2::aes(x = coefficient, y = -log10p, color = p_adjust)) +
      ggplot2::geom_point() +
      ggthemes::theme_clean(base_size = 16) +
      ggplot2::scale_color_gradientn(colors = palette) +
      ggplot2::xlab('Effect Size') +
      ggplot2::ylab('-log10(p)')
    
    # If labels are requested, then add.
    if(label) {
      plot <- plot + ggrepel::geom_text_repel(data = hits[hits$p_adjust < alpha,],
                                              mapping = aes(label = clean_contig_names(Contig_name)),
                                              colour = 'black',
                                              size = 5,
                                              alpha = 0.8,
                                              max.overlaps = Inf)
    }
  }
  
  
  
  return(plot)
}