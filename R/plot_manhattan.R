#' Generate a strain Manhattan plot from a betaGLM object
#'
#' This function creates a Manhattan plot from a `betaGLM` object, with separate facets
#' for the Zero-Inflation and Beta distribution parts of the model.
#'
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate
#' 
#' @param object A `betaGLM` object.
#' @param coef The number of the coefficient from which to generate the plot (default=2).
#' @param taxonomy An optional taxonomy file read using strainspy::read_taxonomy (default=NULL).
#' @param method Character. Multiple testing correction method for p-values (e.g., "holm"). Defaults to "holm".
#' @param alpha Numeric. Significance threshold for adjusted p-values. Defaults to 0.05.
#' @param tax_levels Character vector. Subset of taxonomic levels: c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain") to include in the plot.
#' Defaults to c("Phylum", "Genus", "Species", "Strain"). Provided levels will be reordered by the standard hierarchy and the highest level is used for colouring. If it has >30 categories, 
#' the legend is replaced by a named vector.
#' @param plot If set to false, the function will return a tibble with data used to generate the plot.
#' 
#' @return A `ggplot` object showing the Manhattan plot.
#' 
#' @export
plot_manhattan <- function(object, coef=2, taxonomy=NULL, method = "holm",
                           alpha=0.05, tax_levels = c("Phylum", "Genus", "Species", "Strain"), plot=TRUE) {

  # Validate input
  if (!inherits(object, "betaGLM")) {
    stop("Input must be a betaGLM object.")
  }
  
  # tax_levels and colour_level must be in tax_levels
  tax_levels_ <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
  
  # do some sanity checks
  tax_levels <- intersect(tax_levels_, tax_levels) # drop and reorder
  col_tax_level = tax_levels[1] # choose top level for colouring

  # Create tidy tibbles for plotting
  if (!is.null(taxonomy)) {
    # Add taxonomic informed adjusted p-values
    plot_data <- hadjust(object, coef=coef, taxonomy = taxonomy,
                         method = method, index_range = TRUE)[, c('Level', 'Model', 'Name', 'p_adjust', 'index_min', 'index_max')]

    
    plot_data$Level[!plot_data$Level %in% tax_levels_] <- "Strain"
    plot_data$Level <- factor(plot_data$Level, tax_levels_)
    
    ### issue a warning if the chosen colour has many categories
    legend_ = TRUE
    n_levels = length(table(plot_data$Name[plot_data$Level == col_tax_level]))
    if(n_levels > 30){
      warning(sprintf("There are %d categories of `%s`. Legends >30 are omitted from figure and printed to console.", n_levels, col_tax_level))
      legend_ = FALSE
    }
  
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
      dplyr::filter(Level==col_tax_level) |>
      dplyr::group_by(Name) |>
      dplyr::summarise(
        min_index = min(index_min),
        max_index = max(index_max)
      ) |>
      dplyr::rename(!!col_tax_level := Name)

    
    plot_data <- plot_data |>
      dplyr::cross_join(ranges) |>  # Cross join
      dplyr::filter(index_min >= min_index & index_min <= max_index) |>
      dplyr::select(-min_index, -max_index)  # Remove unnecessary columns

    plot_data <- plot_data |> dplyr::filter(Level %in% tax_levels)

    # Return data if requested
    if (!plot){
      return(plot_data)
    }

    custom_colors = get_colors(n_levels)
    plot_data[[col_tax_level]] <- factor(
      plot_data[[col_tax_level]],
      levels =  unique(plot_data[[col_tax_level]][order(plot_data$index_max)]) # ordering colours in the plot from left to right
    )
   

    gg <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Name, y = log_p_adjust, colour=.data[[col_tax_level]])) +
      ggplot2::geom_segment(ggplot2::aes(y = log_p_adjust, x = index_min - 0.5, xend = index_max + 0.5), size=5) +
      ggplot2::facet_grid(Level ~ Model, scales = "free_x") +
      ggplot2::geom_hline(ggplot2::aes(yintercept = -log10(alpha)), col='red', linetype = "dotted") +
      ggplot2::scale_color_manual(values = custom_colors) +
      ggthemes::theme_clean(base_size = 16) +
      ggplot2::theme(plot.background = ggplot2::element_blank()) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(face = "bold"),
      ) +
      ggplot2::ylab("-log10(p-value)") +
      ggplot2::xlab("Strain")
    
    if(legend_ == FALSE){
        gg <- gg + ggplot2::theme(legend.position = "none")
        categories <- unique(plot_data[[col_tax_level]][order(plot_data$index_max)])
        colors_assigned <- custom_colors[seq_along(categories)]
        names(colors_assigned) <- categories
        cat("Colour Legend:\n")
        print(colors_assigned)
    }

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
