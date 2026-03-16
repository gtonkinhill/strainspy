#' Estimate species richness and coverage of reference strains using ANI
#'
#' This function estimates sample-level **species richness** and 
#' **strain coverage** using ANI-based presence/absence. A strain is 
#' considered present if `ANI > ani_threshold`.
#'
#' Species richness is defined as the number of species with at least
#' one strain exceeding the ANI threshold. This requires taxonomy data, which 
#' can be read in using `read_taxonomy()`.
#'
#' **Experimental** Strain coverage is a database-dependent, normalised score 
#' that increases with the number of reference strains detected within a species. 
#' When only one strain is present for a species, it contributes 1 towards 
#' strain coverage. Asmore strains are detected, the contribution increases 
#' toward 2, reaching 2 when all reference strains are detected. Although this 
#' is designed to prevent inflation for species with many reference strains, 
#' this measures the **coverage of reference strain space in the database**, 
#' not true biological strain richness. Therefore, this value is not returned by
#' default.
#'
#' @param se SummarizedExperiment. A `SummarizedExperiment` object containing the assay data and metadata.
#' @param taxonomy An optional taxonomy file read using `strainspy::read_taxonomy()`.
#' @param ani_threshold Numeric ANI threshold used to define presence (default: 95).
#' @param report_strain_richness Return the estimated strain_richness count. See description (default: FALSE)
#'
#' @return A named list with:
#' \describe{
#'   \item{species_richness}{Integer vector giving the number of species detected per sample.}
#'   \item{strain_coverage}{NA or Numeric vector giving the database-normalised strain coverage per sample.}
#' }
#'
#' @details
#' The strain coverage metric is computed as:
#' \deqn{1 + (k - 1) / N}
#' where \eqn{k} is the number of reference strains detected for a species
#' and \eqn{N} is the total number of reference strains for that species
#' in the database.
#'
#' This approach does **not** resolve true strain multiplicity and may
#' overestimate coverage when a single biological strain matches many
#' highly similar reference genomes. Alternative approaches are being investigated.
#' 
#' @export
estimate_sample_richness = function(se, taxonomy, ani_threshold=95, report_strain_richness = F){
  
  # Get the species level data
  tax_sp = strainspy:::add_tax2tophits(data.frame(Genome_file = se@elementMetadata$Genome_file), taxonomy, columns = "Species")
  
  
  asy = SummarizedExperiment::assay(se)
  species_list <- unique(tax_sp$Species)
  
  present <- Matrix::drop0( (asy >= ani_threshold)*1)  # Essentially a presence/absence representation
  sp_counts <- rowsum(present, group = tax_sp$Species)
  sp_richness = colSums(sp_counts>0)
  
  # Can we do a strain aware approach - measure how much of the strains in the database are covered, but penalise for species with many strains
  if(report_strain_richness){
    N_ref <- table(tax_sp$Species)
    N_ref <- N_ref[rownames(sp_counts)]
    
    N_ref_mat <- matrix(rep(N_ref, ncol(sp_counts)), nrow = nrow(sp_counts), ncol = ncol(sp_counts))
    str_mat <- ifelse(sp_counts == 0, 0, 1 + (sp_counts - 1)/N_ref_mat)
    
    strain_coverage <- colSums(str_mat)
    
  } else {
    strain_coverage <- NA
  }
  
  return(list(species_richness = sp_richness, strain_coverage = strain_coverage))
}

#' Estimate beta diversity from ANI
#'
#' @param se SummarizedExperiment. A `SummarizedExperiment` object containing the assay data and metadata.
#' @param method Distance method: \code{"manhattan"} (default) or \code{"jaccard"}
#' @param distance_only Logical. If TRUE, stop after computing distance matrix
#' @param ani_threshold Numeric. Threshold for presence/absence determination if \code{method == jaccard} (default 95)
#' @param phenotype Any factor variable in `colnames(se@colData)`. Not required if \code{distance_only == FALSE}
#' @param show_plots Logical. If TRUE, shows the plots (default = FALSE)
#' @param return_plots Logical. If TRUE, returns ggplot2 objects (default = FALSE)
#' 
#'
#' @return If distance_only = TRUE, returns distance object. Otherwise a list with:
#' \itemize{
#'   \item \code{distance} — distance matrix
#'   \item \code{nmds} — NMDS object
#'   \item \code{ordination} — NMDS coordinates with group
#'   \item \code{permanova} — PERMANOVA results
#'   \item \code{betadisper} — betadisper object
#'   \item \code{plots} — list of ggplot2 objects (if return_plots = TRUE)
#' }
#'
#' @export
estimate_beta_diversity <- function(
    se,
    method = "manhattan",
    distance_only = FALSE,
    ani_threshold = 95,
    phenotype = NULL,
    show_plots = FALSE,
    return_plots = FALSE
){
  if(distance_only == FALSE & is.null(phenotype)){
    stop("phenotype is required when distance_only = FALSE")
  }
  
  if(!is.null(phenotype)){
    if (!(phenotype %in% colnames(SummarizedExperiment::colData(se)))) {
      stop("phenotype not found in colData(se)")
    }
  }
  
  if(!distance_only | method == "jaccard"){
    if (!requireNamespace("vegan", quietly = TRUE)) {
      stop("The 'vegan' package is required but is not installed.")
    }
  }
  
  # Extract assay
  asy <- t(as.matrix(SummarizedExperiment::assay(se)))
  asy[is.na(asy)] <- 0
  asy <- asy[, apply(asy, 2, sd) > 0.01, drop = FALSE]
  
  # Compute distance
  if(method == "jaccard"){
    asy_bin <- asy >= ani_threshold
    dist_obj <- vegan::vegdist(asy_bin, method = "jaccard")
  } else {
    asy <- 1 - (asy / 100)  # ANI -> genomic distance
    dist_obj <- stats::dist(asy, method = "manhattan") / ncol(asy)
  }
  
  if(distance_only){
    return(dist_obj)
  }
  
  # group
  group <- as.factor(SummarizedExperiment::colData(se)[[phenotype]])
  
  # NMDS
  nmds <- vegan::metaMDS(dist_obj, k = 2)
  ord_df <- as.data.frame(vegan::scores(nmds))
  ord_df$Group <- group
  
  # PERMANOVA
  permanova <- vegan::adonis2(dist_obj ~ group)
  
  # Beta-dispersion
  disp <- vegan::betadisper(dist_obj, group)
  disp_df <- data.frame(
    Sample = rownames(asy),
    Group = group,
    DistanceToCentroid = disp$distances
  )
  
  
  # Get PCoA coordinates from betadisper
  
  plots <- NULL
  
  if(show_plots | return_plots){
    # NMDS ggplot
    nmds_plot <- ggplot2::ggplot(ord_df, ggplot2::aes(NMDS1, NMDS2, color = Group)) +
      ggplot2::geom_point(size = 3) +
      ggplot2::stat_ellipse() +
      ggplot2::theme_bw() +
      ggplot2::labs(title = paste("NMDS (", method, " distance)", sep=""))
    
    # Beta-dispersion boxplot
    betadisp_plot <- ggplot2::ggplot(disp_df, ggplot2::aes(x = Group, y = DistanceToCentroid, fill = Group)) +
      ggplot2::geom_boxplot() +
      ggplot2::theme_bw() +
      ggplot2::labs(title = "Beta-dispersion (distance to group centroid)", y = "Distance to centroid")
    
    plots <- list(
      NMDS = nmds_plot,
      BetaDispersionBox = betadisp_plot
    )
    
  }
  
  if(show_plots){
    cat("Press [Enter] to start visualising plots...")
    # NMDS
    invisible(readline())
    print(plots[[1]])
    
    cat("\nDisplaying plot: NMDS\nPress [Enter] to continue...")
    invisible(readline())
    
    # Beta-dispersion ordination
    print(plots[[2]])
    cat("\nDisplaying plot: Beta-dispersion ordination\nPress [Enter] to continue...")
    invisible(readline())
    
    # Beta-dispersion boxplot
    plot(disp, main = "Beta-dispersion (distance to group centroid)")
    cat("\nDisplaying plot: Beta-dispersion boxplot\nPress [Enter] to continue...")
  }
  
  if(!return_plots) plots = NULL # no one would do this?
  
  
  return(list(
    distance = dist_obj,
    nmds = nmds,
    ordination = ord_df,
    permanova = permanova,
    betadisper = disp,
    plots = plots
  ))
}
