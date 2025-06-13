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
#'
#' }
#' 
#' @export
plot_ani_dist <- function(se, phenotype, contigs, contig_names = NULL, drop_zeros = FALSE, show_points = FALSE, plot = TRUE, plot_type = 'violin'){
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
  
  df_long <- tidyr::pivot_longer(tmp, cols = -phenotype, names_to = "Contig", values_to = "ANI")
  df_long$Contig = factor(df_long$Contig, levels = contig_names)
  
  if(!plot){
    # just return the data.frame
    df_long = as.data.frame(df_long)
    rownames(df_long) = colnames(se)
    return(df_long)
  } else {
    if(drop_zeros){
      df_long$ANI[df_long$ANI == 0] = NA
    }
    df_long$Contig = factor(df_long$Contig)
    phenotype_ <- rlang::sym(phenotype)
    
    p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Contig, y = ANI, fill = !!phenotype_))
    
    if(plot_type == 'violin'){
      p <- p + ggplot2::geom_violin(trim = T, alpha = 0.2, drop = FALSE)
    } else {
      p <- p + ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.5)  # semi-transparent boxplot
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