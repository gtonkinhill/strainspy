#' Generate a violin or box plot for a subset of contigs for a categorical phenotype.
#' @importFrom SummarizedExperiment assay
#' @importFrom tidyr pivot_longer
#' @importFrom rlang sym
#' @importFrom ggplot2 ggplot aes geom_violin theme_minimal labs scale_fill_brewer theme element_text
#' @importFrom ggforce geom_sina
#'
#' @param se  SummarizedExperiment. A `SummarizedExperiment` object containing the assay data and metadata.
#' @param phenotype One variable in `colnames(se@colData)`, except `Sample_File`. The selected variable must be categorical
#' @param contigs A vector of contigs to generate violin plots
#' @param contig_names Optional character vector to replace long contig names in the plot.
#' Must match the length and order of `contigs`. If NULL (default), shortening will be attempted automatically.
#' @param drop_zeros Bool. If TRUE, ANI or Abundance 0 values will not be included in the violin plot. Default FALSE.
#' @param show_points Bool. If TRUE, points will be overlaid with jitter using ggforce::geom_sina(). Default FALSE.
#' @param plot If set to false, the function will return a tibble with data used to generate the plot.
#' @param plot_type Specify whether `violin` or `box` plot is required. Default `violin`.
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
#'
#'
#' }
#'
#' @export
plot_ani_dist <- function(se, phenotype, contigs, contig_names = NULL, drop_zeros = FALSE, 
                          show_points = FALSE, plot = TRUE, plot_type = 'violin'){
  if(length(phenotype) != 1) {
    stop("Only one phenotype can be plotted at a time")
  }
  
  if(! (phenotype %in% colnames(SummarizedExperiment::colData(se))) ){
    stop("Phenotype not found in colnames(se@colData)")
  }
  
  if( ! is.factor(se@colData[[phenotype]])) {
    stop("Phenotype must be a factor")
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
    plot_type = 'violin'
  }
  
  tmp = as.data.frame(t(as.matrix(SummarizedExperiment::assay(se[contigs, ]))))
  tmp = cbind(se@colData[[phenotype]], tmp)
  
  if(!is.null(contig_names)){
    colnames(tmp) = c(phenotype, contig_names)
  } else {
    contig_names = clean_contig_names(colnames(tmp)[-1])
    colnames(tmp) = c(phenotype, contig_names)
  }
  
  tmp <- tibble::rownames_to_column(tmp, var = "sample")
  df_long <- tidyr::pivot_longer(tmp, cols = -c(sample, phenotype), names_to = "Contig", values_to = "ANI")
  df_long$Contig = factor(df_long$Contig, levels = contig_names)
  
  if(!plot){
    # just return the data.frame
    return(df_long)
  } else {
    if(drop_zeros){
      df_long$ANI[df_long$ANI == 0] = NA
    }
    df_long$Contig = factor(df_long$Contig)
    phenotype_ <- rlang::sym(phenotype)
    
    p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Contig, y = ANI, group = interaction(Contig, !!phenotype_)))
    
    if(plot_type == 'violin'){
      p <- p + ggplot2::geom_violin(ggplot2::aes(fill = !!phenotype_),
                                    trim = T, alpha = 0.2, drop = FALSE)
    } else {
      p <- p + ggplot2::geom_boxplot(ggplot2::aes(fill = !!phenotype_),
                                     outlier.shape = NA, alpha = 0.5)  # semi-transparent boxplot
    }
    
    
    p <- p +
      ggplot2::labs(x = "Contigs", y = "ANI") +
      ggplot2::scale_fill_brewer(palette = "Set2") +
      ggplot2::scale_color_brewer(palette = "Set2") +
      ggplot2::theme_minimal(base_size = 16) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    
    if (show_points) {
      p <- p + ggforce::geom_sina(
        ggplot2::aes(color = !!phenotype_),
        alpha = 1, na.rm = TRUE, method = "count", seed = 1988
      )
    }
    print(p)
    
  }
  
  
}

quick_histo <- function(se, phenotype, contig, drop_zeros = F) {
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
