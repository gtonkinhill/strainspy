#' Generate a Violin plot for a subset of contigs for a given categorical phenotype 
#' @importFrom SummarizedExperiment assay
#' @importFrom tidyr pivot_longer
#' @importFrom rlang sym
#' @importFrom ggplot2 ggplot aes geom_violin theme_minimal labs scale_fill_brewer theme element_text
#' 
#' @param se  SummarizedExperiment. A `SummarizedExperiment` object containing the assay data and metadata.
#' @param phenotype One variable in `colnames(se@colData)`, except `Sample_File`. The selected variable must be categorical
#' @param contigs A vector of contigs to generate violin plots
#' @param contig_names A vector of names for the contigs to use in the plot, useful when the original contig names are too long. Must be the same length and order as `contig_names`. Default NULL (shorten automatically).
#' @param drop_zeros Bool. If TRUE, ANI or Abundance 0 values will not be included in the violin plot. Default FALSE.
#' @param plot If set to false, the function will return a tibble with data used to generate the plot.
#' 
#' @export
plot_violin <- function(se, phenotype, contigs, contig_names = NULL, drop_zeros = FALSE, plot = TRUE){
  if(length(phenotype) != 1) {
    stop("Only one phenotype can be plotted at a time")
  }
  
  if(! (phenotype %in% colnames(se@colData)) ){
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
  
  tmp = as.data.frame(t(as.matrix(SummarizedExperiment::assay(se[contigs, ]))))
  tmp = cbind(se@colData[[phenotype]], tmp)
  
  if(!is.null(contig_names)){
    colnames(tmp) = c(phenotype, contig_names)
  } else {
    colnames(tmp) = c(phenotype, clean_contig_names(colnames(tmp)[-1]))
  }
  
  df_long <- tidyr::pivot_longer(tmp, cols = -phenotype, names_to = "Contig", values_to = "ANI")
  
  if(drop_zeros){
    df_long$ANI[df_long$ANI == 0] = NA
  }
  
  phenotype <- rlang::sym(phenotype)
  
  if(!plot){
    return(df_long)
  } else {
    ggplot(df_long, aes(x = Contig, y = ANI, fill = !!phenotype)) +
      geom_violin(trim = FALSE, alpha = 0.6) +
      theme_minimal() +
      labs(x = "Contigs", y = "ANI") +
      scale_fill_brewer(palette = "Set2") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  
}