#' Generate a strain Manhattan plot from a betaGLM object
#'
#' This function creates a Manhattan plot from a `betaGLM` object, with separate facets
#' for the Zero-Inflation and Beta distribution parts of the model.
#'
#' @param object A `betaGLM` object.
#' @param coef The number of the coefficient from which to generate the plot (default=2).
#' @param plot If set to false, the function will return a tibble with data used to generate the plot.
#' @return A `ggplot` object showing the Manhattan plot.
#' @export
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate
plot_manhattan <- function(object, coef=2, taxonomy=NULL, method = "HMP",
                           alpha=0.05, levels = c("Phylum", "Genus", "Species", "Strain"),
                           colour_level = "Phylum", plot=TRUE) {

  # Validate input
  if (!inherits(object, "betaGLM")) {
    stop("Input must be a betaGLM object.")
  }

  # Create tidy tibbles for plotting
  if (!is.null(taxonomy)) {
    # Add taxonomic informed adjusted p-values
    plot_data <- hadjust(object, coef=coef, taxonomy = taxonomy,
                         method = method, index_range = TRUE)[, c('Level', 'Model', 'Name', 'p_adjust', 'index_min', 'index_max')]

    tax_levels <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
    plot_data$Level[!plot_data$Level %in% tax_levels] <- "Strain"
    plot_data$Level <- factor(plot_data$Level, levels=tax_levels)
  } else {
    # Add sequence level adjusted p-values
    hits <- top_hits(object, coef=coef, method = method, alpha = 1)

    plot_data <- hits[,c(colnames(object@row_data), "p_adjust")] |>
      tibble::add_column(Model="Beta", .before=1)

    if ('zi_p_adjust' %in% colnames(hits)) {
      plot_data <- rbind(plot_data,
                         hits[,c(colnames(object@row_data), "zi_p_adjust")] |>
                           dplyr::rename(p_adjust = "zi_p_adjust") |>
                           tibble::add_column(Model="Zero-Inflated", .before=1))
    }

    colnames(plot_data)[[2]] <- "Name"
    plot_data <- plot_data |> tibble::add_column(Level="Sequence", .before=1)

    plot_data$index <- 1:nrow(plot_data)
  }

  plot_data$log_p_adjust <- -log10(plot_data$p_adjust)


  # Create the Manhattan plot
  if (!is.null(taxonomy)){
    # add colours for chosen level
    ranges <- plot_data |>
      dplyr::filter(Level==colour_level) |>
      dplyr::group_by(Name) |>
      dplyr::summarise(
        min_index = min(index_min),
        max_index = max(index_max)
      ) |>
      dplyr::rename(Phylum = Name)

    plot_data <- plot_data |>
      dplyr::cross_join(ranges) |>  # Cross join
      dplyr::filter(index_min >= min_index & index_min <= max_index) |>
      dplyr::select(-min_index, -max_index)  # Remove unnecessary columns

    plot_data <- plot_data |> dplyr::filter(Level %in% levels)

    # Return data if requested
    if (!plot){
      return(plot_data)
    }

    # Define your custom color palette
    custom_colors <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#b15928',
                       '#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5')

    custom_colors <- custom_colors[1:length(unique(plot_data$Phylum))]

    gg <- ggplot2::ggplot(plot_data, aes(x = Name, y = log_p_adjust, colour=Phylum)) +
      ggplot2::geom_segment(aes(y = log_p_adjust, x = index_min - 0.5, xend = index_max + 0.5), size=5) +
      ggplot2::facet_grid(Level ~ Model, scales = "free_x") +
      geom_hline(aes(yintercept = -log10(alpha)), col='red', linetype = "dotted") +
      scale_color_manual(values = custom_colors) +
      ggthemes::theme_clean(base_size = 16) +
      ggplot2::theme(plot.background = element_blank()) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(face = "bold"),
      ) +
      ggplot2::ylab("-log10(p-value)") +
      ggplot2::xlab("Strain")

  } else {

    # Return data if requested
    if (!plot){
      return(plot_data)
    }

    gg <- ggplot2::ggplot(plot_data, aes(x = Name, y = log_p_adjust)) +
      ggplot2::geom_point() +
      ggplot2::facet_wrap( ~ Model) +
      geom_hline(aes(yintercept = -log10(alpha)), col='red', linetype = "dotted") +
      ggthemes::theme_clean(base_size = 16) +
      ggplot2::theme(plot.background = element_blank()) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(face = "bold"),
      ) +
      ggplot2::ylab("-log10(p-value)") +
      ggplot2::xlab("Strain")

  }

  return(gg)
}

# Example usage
# Assuming `beta_glm_obj` is a valid betaGLM object
# plot_manhattan(beta_glm_obj)
