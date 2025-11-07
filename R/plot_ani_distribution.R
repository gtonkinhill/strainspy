#' Generate a violin or box plot for a subset of contigs for a categorical phenotype.
#' @importFrom SummarizedExperiment assay
#' @importFrom tidyr pivot_longer
#' @importFrom rlang sym
#' @importFrom ggplot2 ggplot aes geom_violin theme_minimal labs scale_fill_brewer theme element_text
#' @importFrom ggforce geom_sina
#' @importFrom patchwork plot_layout
#' 
#' 
#' @param se  SummarizedExperiment. A `SummarizedExperiment` object containing the assay data and metadata.
#' @param phenotype Any factor variable in `colnames(se@colData)`.
#' @param contigs A vector of contigs to generate distribution plots
#' @param facet_phenotype Any factor variable in `colnames(se@colData)`, except `phenotype`. If provided, this will be used to facet the plot. See `?ggplot2::facet_wrap()`. Default NULL.
#' @param contig_names Optional character vector to replace long contig names in the plot.
#' Must match the length and order of `contigs`. If NULL (default), shortening will be attempted automatically.
#' @param show_zero_plot Bool. If TRUE, distribution of the proportion of zeros will also be plotted. Default TRUE
#' @param show_points Bool. If TRUE, points will be overlaid with jitter using ggforce::geom_sina(). Default TRUE
#' @param plot If set to false, the function will return a tibble with data used to generate the plot.
#' @param plot_type Specify whether `box` or `violin` plot is required. Default `box`.
#'
#' @examples
#' \dontrun{
#' library(strainspy)
#'
#' example_meta_path <- system.file("extdata", "example_metadata.csv.gz", package = "strainspy")
#' example_meta <- readr::read_csv(example_meta_path)
#' # Change variable from character to factor
#' example_meta$Case_status = factor(example_meta$Case_status)
#'
#' example_path <- system.file("extdata", "example_sylph_profile.tsv.gz", package = "strainspy")
#' se <- read_sylph(example_path, example_meta)
#'
#' plot_ani_dist(se, phenotype = "Case_status", contigs = rownames(se)[1:5])
#' }
#'
#' @export
plot_ani_dist <- function(se, phenotype, contigs, facet_phenotype = NULL, contig_names = NULL, show_zero_plot = TRUE, 
                          show_points = TRUE, plot = TRUE, plot_type = 'box'){
  if(length(phenotype) != 1) {
    stop("Only one phenotype can be plotted at a time")
  }
  
  if(! (phenotype %in% colnames(SummarizedExperiment::colData(se))) ){
    stop("phenotype not found in colnames(se@colData)")
  }
  
  if( ! is.factor(se@colData[[phenotype]])) {
    stop("phenotype must be a factor")
  }
  
  if(!is.null(facet_phenotype)) {
    facet = T
    
    if(phenotype == facet_phenotype){
      stop("phenotype and facet_phenotype cannnot be the same")
    }
    
    if(length(facet_phenotype) != 1) {
      stop("Only one facet_phenotype can be used at a time")
    } 
    
    if(! (facet_phenotype %in% colnames(SummarizedExperiment::colData(se))) ){
      stop("facet_phenotype not found in colnames(se@colData)")
    }
    
    if( ! is.factor(se@colData[[facet_phenotype]])) {
      stop("facet_phenotype must be a factor")
    }
  } else {
    facet = F
  }
  
  chk_contigs = contigs %in% rownames(se)
  if( ! all(chk_contigs)  ){
    cat(paste(contigs[!chk_contigs], '\n', sep = ""))
    stop("Above contigs are missing from rownames(se)")
  }
  
  if (!is.null(contig_names) && length(contig_names) != length(contigs)) {
    stop("`contig_names` must be the same length as `contigs`")
  }
  
  if (!plot_type %in% c('violin', 'box')) {
    plot_type = 'box'
  }
  
  # Lazy coding, fix later - shouldn't impact speed
  if(!facet){
    tmp = as.data.frame(t(as.matrix(SummarizedExperiment::assay(se[contigs, ]))))
    tmp = cbind(se@colData[[phenotype]], tmp)
    
    if(!is.null(contig_names)){
      colnames(tmp) = c(phenotype, contig_names)
    } else {
      contig_names = strainspy:::clean_contig_names(colnames(tmp)[-1])
      colnames(tmp) = c(phenotype, contig_names)
    }
  } else {
    tmp = as.data.frame(t(as.matrix(SummarizedExperiment::assay(se[contigs, ]))))
    tmp = cbind(se@colData[[phenotype]], se@colData[[facet_phenotype]], tmp)
    
    if(!is.null(contig_names)){
      colnames(tmp) = c(phenotype, facet_phenotype, contig_names)
    } else {
      contig_names = strainspy:::clean_contig_names(colnames(tmp)[-c(1,2)])
      colnames(tmp) = c(phenotype, facet_phenotype, contig_names)
    }
  }
  
  
  tmp <- tibble::rownames_to_column(tmp, var = "sample")
  
  
  if(!facet){
    df_long <- tidyr::pivot_longer(
      tmp,
      cols = -all_of(c("sample", phenotype)),
      names_to = "Contig",
      values_to = "ANI"
    )
  } else {
    df_long <- tidyr::pivot_longer(
      tmp,
      cols = -all_of(c("sample", phenotype, facet_phenotype)),
      names_to = "Contig",
      values_to = "ANI"
    )
  }
  
  df_long$Contig = factor(df_long$Contig, levels = contig_names)
  
  # Handle NAs first
  cmpcases = complete.cases(df_long)
  
  if(any(cmpcases == FALSE)){
    df_long = df_long[which(cmpcases), ]
    
    if(nrow(df_long) == 0) {
      stop("All entries contain incomplete entries, unable to plot")
    } else {
      warning("There are incomplete entries in: ", length(which(cmpcases == FALSE)), " samples, these will be dropped")  
    }
    
  }
  
  
  # Get the zeros out
  if(!facet) {
    zero_props <- df_long %>%
      mutate(!!sym(phenotype) := as.factor(.data[[phenotype]])) %>%
      group_by(Contig, !!sym(phenotype)) %>%
      summarise(
        zero_count = sum(ANI == 0),
        total = n(),
        zero_prop = zero_count / total,
        .groups = "drop"
      )
  } else {
    zero_props <- df_long %>%
      mutate(!!sym(phenotype) := as.factor(.data[[phenotype]])) %>%
      group_by(Contig, !!sym(phenotype), !!sym(facet_phenotype)) %>%
      summarise(
        zero_count = sum(ANI == 0),
        total = n(),
        zero_prop = zero_count / total,
        .groups = "drop"
      )
  }
  
  if(!plot){
    # just return the data.frame
    return(df_long)
  } else {
    df_long$ANI[df_long$ANI == 0] = NA
    
    df_long$Contig = factor(df_long$Contig)
    phenotype_ <- rlang::sym(phenotype)
    
    phen_levels <- levels(factor(df_long[[phenotype]]))
    phen_cols <- setNames(strainspy:::get_colors(length(phen_levels)), phen_levels) # can maintain consistency between the two plots
    
    p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Contig, y = ANI, group = interaction(Contig, !!phenotype_)))
    
    if(plot_type == 'violin'){
      p <- p + ggplot2::geom_violin(ggplot2::aes(fill = !!phenotype_),
                                    trim = T, alpha = 0.2, drop = FALSE)
    } else {
      p <- p + ggplot2::geom_boxplot(ggplot2::aes(fill = !!phenotype_),
                                     outlier.shape = NA, alpha = 0.5)  # semi-transparent boxplot
    }
    
    
    p <- p +
      ggplot2::labs(y = "ANI", x = NULL) +
      ggplot2::theme_minimal(base_size = 16)
    
    if (show_points) {
      p <- p + ggforce::geom_sina(
        ggplot2::aes(color = !!phenotype_),
        alpha = 1, na.rm = TRUE, method = "count", seed = 1988
      )
    }
    
    p <- p + ggplot2::scale_fill_manual(values = phen_cols, name = phenotype) +
      ggplot2::scale_color_manual(values = phen_cols, name = phenotype)
    
    
    if(facet){
      p <- p + ggplot2::facet_wrap(as.formula(paste("~", facet_phenotype))) +
        ggplot2::ggtitle(paste("Faceted by", facet_phenotype))
    }
    
    if(show_zero_plot){
      p_p = list()
      p_p[[1]] <- p  +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(),       
                       axis.ticks.x = ggplot2::element_blank(),      
                       axis.line.x = ggplot2::element_blank())       # remove x-axis
      
      p2 <- ggplot2::ggplot(zero_props, ggplot2::aes(x = Contig, y = zero_prop, fill = !!sym(phenotype))) +
        ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(width = 0.4), width = 0.35) +
        ggplot2::theme_minimal(base_size = 16) +
        ggplot2::labs(x = "Contig", y = "Absence proportion", fill = phenotype) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          panel.grid.major.x = ggplot2::element_blank()
        ) + 
        ggplot2:: scale_fill_manual(values = phen_cols, name = phenotype)
      
      if(facet){
        p2 <- p2 +  ggplot2::facet_wrap(as.formula(paste("~", facet_phenotype)))
      }
      
      p_p[[2]] <- p2
      
      return(p_p[[1]] / p_p[[2]] + patchwork::plot_layout(heights = c(3, 1)))
    } else {
      return(p)
    }
    
    
    
    
  }
}

#' Generate a histogram for a single contig for a categorical phenotype.
#' @importFrom SummarizedExperiment assay
#' @importFrom tidyr pivot_longer
#' @importFrom rlang sym
#' @importFrom ggplot2 ggplot aes geom_violin theme_minimal labs scale_fill_brewer theme element_text
#' @importFrom ggforce geom_sina
#'
#' @param se  SummarizedExperiment. A `SummarizedExperiment` object containing the assay data and metadata.
#' @param phenotype One variable in `colnames(se@colData)`, except `Sample_File`. The selected variable must be categorical
#' @param contig Contig name to visualise
#' @param drop_zeros Bool. If TRUE, ANI or Abundance 0 values will not be included in the violin plot. Default FALSE
#'
#' @examples
#' \dontrun{
#' library(strainspy)
#'
#' example_meta_path <- system.file("extdata", "example_metadata.csv.gz", package = "strainspy")
#' example_meta <- readr::read_csv(example_meta_path)
#' # Change variable from character to factor
#' example_meta$Case_status = factor(example_meta$Case_status)
#'
#' example_path <- system.file("extdata", "example_sylph_profile.tsv.gz", package = "strainspy")
#' se <- read_sylph(example_path, example_meta)
#'
#' plot_histogram(se, phenotype = "Case_status", contigs = rownames(se)[1])
#'
#' }
#'
#' @export
plot_histogram <- function(se, phenotype, contig, drop_zeros = F) {
  stopifnot(is(se, "SummarizedExperiment"))
  stopifnot(is.character(phenotype), length(phenotype) == 1)
  
  dat <- colData(se) |> as.data.frame()
  
  idx = which(rownames(se) == contig)
  if (length(idx) != 1) {
    stop(sprintf("Contig '%s' not found (or found multiple times) in rowNames(se).", contig))
  }
  
  dat$Value_orig <- assay(se)[idx, ]
  
  
  if (drop_zeros) {
    dat <- dplyr::filter(dat, Value_orig != 0)
  }
  
  summ <- dat |>
    dplyr::group_by(.data[[phenotype]]) |>
    dplyr::summarise(
      Min    = min(Value_orig),
      Q1     = quantile(Value_orig, 0.25),
      Median = median(Value_orig),
      Mean   = mean(Value_orig),
      Q3     = quantile(Value_orig, 0.75),
      Max    = max(Value_orig),
      .groups = "drop"
    ) |>
    tidyr::pivot_longer(
      -dplyr::all_of(phenotype),
      names_to = "stat",
      values_to = "value"
    )
  
  p <- ggplot2::ggplot(dat,
                       ggplot2::aes(x = Value_orig, fill = .data[[phenotype]])) +
    ggplot2::geom_histogram(position = "identity", alpha = 0.4, bins = 250) +
    ggplot2::geom_vline(data = summ,
                        ggplot2::aes(xintercept = value, color = stat),
                        linetype = "dashed", size = 0.6, show.legend = TRUE) +
    ggplot2::facet_wrap(stats::as.formula(paste("~", phenotype)), ncol = 1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Value_orig", y = "Count", color = "Statistic") +
    ggplot2::ggtitle(clean_contig_names(contig)) +
    ggplot2::scale_color_manual(values = c(
      Min    = "black",
      Q1     = "blue",
      Median = "red",
      Mean   = "darkgreen",
      Q3     = "blue",
      Max    = "black"
    )) +
    ggplot2::theme(text = ggplot2::element_text(size = 16))
  
  attr(p, "summary_table") <- summ
  return(p)
}
