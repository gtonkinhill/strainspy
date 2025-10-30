#' Generate a strain Manhattan plot from a strainspy_fit object
#'
#' This function can generate multiple Manhattan plots from a `strainspy_fit` object. 
#' 
#' The generated plot will be more informative if taxonomic information is 
#' provided. Two types of plots can be generated with taxonomic information.
#' 
#' 1. Set `aggregate_by_taxa = TRUE` to aggregate p-values across taxonomic levels 
#' specified by `tax_levels`. This generates a stacked segment Manhattan plot. 
#' In this plot, -log10 p-values adjusted using `method` will be plotted in 
#' taxonomic ordering. Significance threshold will be `alpha`.
#' 
#' 2. Set `aggregate_by_taxa = FALSE` to generate a traditional Manhattan plot for 
#' all strains coloured by phylum. Unadjusted p-values will be plotted and 
#' Bonferroni significance thresholds (0.05, 0.01) will be shown. Strains will 
#' be in taxonomic ordering and a tree showing the taxonomic relationship between 
#' strains will be plotted below the Manhattan plot.
#' 
#' If taxonomic information is not provided (`taxonomy = NULL`), a plot with 
#' -log10 p-values adjusted using `method` will be generated. Strains will be 
#' plotted in alphabetical order. Significance threshold will be `alpha`.
#' 
#' If Zero-Inflated Beta regression was used to model `object`, separate 
#' facets will be used for zero-inflated and beta components of the model. 
#'
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate
#' 
#' @param object A `strainspy_fit` object.
#' @param coef The number of the coefficient from which to generate the plot (default=2).
#' @param taxonomy An optional taxonomy file read using strainspy::read_taxonomy (default=NULL).
#' @param aggregate_by_taxa Logical. If TRUE, aggregate p-values across `tax_levels` and visualise as a stacked segment plot. Requires `taxonomy`. If NULL (default), the function automatically sets it to TRUE when `taxonomy` is provided.
#' @param method Character. Multiple testing correction method for p-values (e.g., "holm"). Defaults to "holm". Only applicable if `aggregate_by_taxa = TRUE`.
#' @param alpha Numeric. Significance threshold for adjusted p-values. Defaults to 0.05.
#' @param tax_levels Character vector. Subset of taxonomic levels: c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain") to include in the plot.
#' Defaults to c("Phylum", "Genus", "Species", "Strain"). Provided levels will be reordered by the standard hierarchy and the highest level is used for colouring. If it has >30 categories, 
#' the legend is replaced by a named vector. Only applicable with `aggregate_by_taxa = TRUE`.
#' @param plot If set to false, the function will return a tibble with data used to generate the plot. NOTE: Will be disregarded if taxonomy is provided with the option `aggregate_by_taxa = FALSE`.
#' 
#' @return A `ggplot` object showing the Manhattan plot if `plot = TRUE`, else a `data.frame` with plot data.
#' 
#' @export
plot_manhattan <- function(object, coef=2, taxonomy=NULL, aggregate_by_taxa = NULL, method = "holm",
                           alpha=0.05, tax_levels = c("Phylum", "Genus", "Species", "Strain"), plot=TRUE) {
  
  # sanity checks
  # Validate input
  if (!inherits(object, "strainspy_fit")) {
    stop("Input must be a strainspy_fit object.")
  }
  
  # Handle default for aggregate_by_taxa
  if (is.null(aggregate_by_taxa)) {
    aggregate_by_taxa <- !is.null(taxonomy)
  }
  
  # Consistency check
  if (isTRUE(aggregate_by_taxa) && is.null(taxonomy)) {
    stop("Cannot aggregate by taxa without providing taxonomy.")
  }
  
  # tax_levels and colour_level must be in tax_levels
  tax_levels_ <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
  
  # do some sanity checks
  tax_levels <- intersect(tax_levels_, tax_levels) # drop and reorder
  
  if(!"Strain" %in% tax_levels) tax_levels = c(tax_levels, "Strain")
  
  col_tax_level = tax_levels[1] # choose top level for colouring
  
  
  if (!is.null(taxonomy)) {
    
    if(isFALSE(aggregate_by_taxa)){ # requesting the strain level MHP with tree
      return(plot_manhattan_tree(object, taxonomy, coef = coef))
    }
    
    
    # Create tidy tibbles for plotting
    # Add taxonomic informed adjusted p-values
    plot_data <- hadjust(object, coef=coef, taxonomy = taxonomy,
                         method = method, index_range = TRUE)[, c('Level', 'Model', 'Name', 'p_adjust', 'index_min', 'index_max')]
    
    
    plot_data$Level[!plot_data$Level %in% tax_levels_] <- "Strain" # Rename contig_names to strain
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
    
    custom_colors = strainspy:::colour_by_tax(genomes = plot_data$Name[which(plot_data$Level == "Strain")], taxonomy = taxonomy,
                                  tax_levels = tax_levels)
    plot_data[[col_tax_level]] <- factor(
      plot_data[[col_tax_level]],
      levels = names(custom_colors)
      # levels =  unique(plot_data[[col_tax_level]][order(plot_data$index_max)]) # ordering colours in the plot from left to right
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
        panel.spacing = unit(1.5, "lines"),
        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5)
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

#' Generate a `conventional` Manhattan plot for the strain all strains. 
#'
#' @importFrom ggtree ggtree %<+% theme_tree2 geom_tippoint
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual geom_hline labs theme element_text element_blank
#' @importFrom ggthemes theme_clean
#' @importFrom dplyr filter arrange pull
#' @importFrom ape as.phylo
#' @importFrom patchwork plot_layout
#'
#' @param object A `strainspy_fit` object.
#' @param taxonomy Taxonomy file read using strainspy::read_taxonomy (default=NULL).
#' @param coef The number of the coefficient from which to generate the plot (default=2).
#'
#' @return A `ggplot` object showing the Manhattan plot combined with the taxonomic tree
plot_manhattan_tree <- function(object, taxonomy, coef = 2) {
  
  if (!inherits(object, "strainspy_fit")) {
    stop("Input must be a strainspy_fit object.")
  }
  
  # tax_levels and colour_level must be in tax_levels
  tax_levels_ <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
  # tax_levels = c("Domain", "Phylum", "Genus", "Species", "Strain")
  # # do some sanity checks
  # tax_levels <- intersect(tax_levels_, tax_levels) # drop and reorder
  # col_tax_level = tax_levels[1] # choose top level for colouring
  
  # Create tidy tibbles for plotting
  # Add taxonomic informed adjusted p-values
  plot_data_ <- hadjust(object, coef=coef, taxonomy = taxonomy,
                        method = "holm", index_range = TRUE)[, c('Level', 'Model', 'Name', 'p_adjust', 'index_min', 'index_max')]
  
  
  # number of models we need to plot for (only applies to zib/beta vs others for now)
  models = unique(plot_data_$Model)
  p_tree = list()
  p_mh = list()
  
  for(m in models){
    plot_data = plot_data_[which(plot_data_$Model == m), ]
    
    plot_data$Level[!plot_data$Level %in% tax_levels_] <- "Strain" # Rename contig_names to strain
    plot_data$Level <- factor(plot_data$Level, tax_levels_)
    
    if ("Strain" %in% plot_data$Level) {
      var <- "Strain"
    } else {
      if(!"Contig_name" %in% plot_data$Level){
        stop("Strain or Contig_name required to generate plot.")
      }
      var <- "Contig_name"
    }
    
    idx_ord <- match(plot_data$Name[plot_data$Level == var], object@row_data@listData$Genome_file)
    if(m == "Zero-Inflated"){
      p_vals_orig <- object@zi_p_values@listData[[coef]][idx_ord]
    } else {
      p_vals_orig <- object@p_values@listData[[coef]][idx_ord]
    }
    
    
    dat_mhp <- data.frame(
      Genome = object@row_data@listData$Genome_file[idx_ord],
      p_vals_orig = -log10(p_vals_orig)
    )
    
    ## let's ensure alphabetic ordering for Genome here
    # The tree can sometimes differ based on genome order
    dat_mhp = dat_mhp[order(dat_mhp$Genome),]; rownames(dat_mhp) = NULL
    
    tax_mhp <- taxonomy[match(dat_mhp$Genome, taxonomy$Genome), ]
    tax_mhp <- tax_mhp[order(tax_mhp$Domain), ]
    
    tax_mhp <- data.frame(
      Genome = factor(tax_mhp$Genome),
      Domain = factor(tax_mhp$Domain),
      Phylum = factor(tax_mhp$Phylum),
      Class = factor(tax_mhp$Class),
      Order = factor(tax_mhp$Order),
      Family = factor(tax_mhp$Family),
      Genus = factor(tax_mhp$Genus)
    )
    
    tree <- ape::as.phylo.formula(~ Domain/Phylum/Class/Order/Family/Genus/Genome, data = tax_mhp)
    
    # Pre-plot to extract tip order
    p_tmp <- ggtree::ggtree(tree)
    tip_order <- p_tmp$data |>
      dplyr::filter(isTip) |>
      dplyr::arrange(y) |>
      dplyr::pull(label)
    
    tax_mhp <- tax_mhp[match(tip_order, tax_mhp$Genome), ]; rownames(tax_mhp) <- NULL
    dat_mhp <- dat_mhp[match(tip_order, dat_mhp$Genome), ]; rownames(dat_mhp) <- NULL
    
    phylum_order <- unique(tax_mhp$Phylum)
    tax_mhp$Phylum <- factor(tax_mhp$Phylum, levels = phylum_order)
    dat_mhp$Phylum <- tax_mhp$Phylum
    dat_mhp$Genome <- factor(dat_mhp$Genome, levels = dat_mhp$Genome)
    
    color_palette <- setNames(
      get_colors(length(phylum_order)),
      phylum_order
    )
    
    p_tree[[m]] <- ggtree(tree, layout = "rectangular") %<+% tax_mhp +
      geom_tippoint(aes(color = Phylum), size = 2) +
      ggplot2::coord_flip() +
      ggplot2::scale_color_manual(values = color_palette) +
      theme_tree2() +
      ggplot2::theme(legend.position = "bottom", text = ggplot2::element_text(size = 16))
    
    p_mh[[m]] <- ggplot2::ggplot(dat_mhp, aes(x = Genome, y = p_vals_orig, color = Phylum)) +
      ggplot2::geom_point(size = 2) +
      ggplot2::scale_color_manual(values = color_palette) +
      ggplot2::geom_hline(aes(yintercept = -log10(0.05 / nrow(dat_mhp))), col = 'black', linetype = "dashed") +
      ggplot2::geom_hline(aes(yintercept = -log10(0.01 / nrow(dat_mhp))), col = 'red', linetype = "dashed") +
      ggthemes::theme_clean(base_size = 16) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        legend.position = "none"
      ) +
      ggplot2::labs(x = NULL, y = "-log10(p)") +
      ggplot2::ggtitle(m)
    
  }
  
  # generate the final plot
  # Ensure both lists are same length
  stopifnot(length(p_mh) == length(p_tree))
  
  # This one is fancy, but there is no point in plotting the tree twice!
  # # Combine each Manhattan/tree pair vertically
  # combined_pairs <- Map(function(mh, tree, title) {
  #   mh <- mh + ggplot2::ggtitle(title)
  #   mh / tree + patchwork::plot_layout(heights = c(4, 1))
  # }, p_mh, p_tree, models)
  
  # Now lay them out side by side
  # return(patchwork::wrap_plots(combined_pairs, nrow = 1))
  
  if(length(models) == 2) {
    return(p_mh[[1]] / p_mh[[2]] / p_tree[[1]] + patchwork::plot_layout(heights = c(4, 4, 1)))
  } else {
    return(p_mh[[1]] / p_tree[[1]] + patchwork::plot_layout(heights = c(4, 1)))
  }
  
}
