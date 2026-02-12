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
#' Strain coverage is a database-dependent, normalised score that increases
#' with the number of reference strains detected within a species. When only one 
#' strain is present for a species, it contributes 1 towards strain coverage. As 
#' more strains are detected, the contribution increases toward 2, reaching 2 
#' when all reference strains are detected. This is designed to prevent 
#' inflation for species with many reference strains.
#'
#' This measures the **coverage of reference strain space in the database**, 
#' not true biological strain richness.
#'
#' @param se SummarizedExperiment. A `SummarizedExperiment` object containing the assay data and metadata.
#' @param taxonomy An optional taxonomy file read using `strainspy::read_taxonomy()`.
#' @param ani_threshold Numeric ANI threshold used to define presence (default: 95).
#'
#' @return A named list with:
#' \describe{
#'   \item{species_richness}{Integer vector giving the number of species detected per sample.}
#'   \item{strain_coverage}{Numeric vector giving the database-normalised strain coverage per sample.}
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
#' highly similar reference genomes.
estimate_sample_richness = function(se, taxonomy,ani_threshold=95){
  
  # Get the species level data
  tax_sp = strainspy:::add_tax2tophits(data.frame(Genome_file = se@elementMetadata$Genome_file), taxonomy, columns = "Species")
  
  
  asy = SummarizedExperiment::assay(se)
  species_list <- unique(tax_sp$Species)
  
  present <- Matrix::drop0( (asy >= ani_threshold)*1)  # Essentially a presence/absence representation
  sp_counts <- rowsum(present, group = tax_sp$Species)
  sp_richness = colSums(sp_counts>0)
  
  # Can we do a strain aware approach - measure how much of the strains in the database are covered, but penalise for species with many strains
  N_ref <- table(tax_sp$Species)
  N_ref <- N_ref[rownames(sp_counts)]

  N_ref_mat <- matrix(rep(N_ref, ncol(sp_counts)), nrow = nrow(sp_counts), ncol = ncol(sp_counts))
  str_mat <- ifelse(sp_counts == 0, 0, 1 + (sp_counts - 1)/N_ref_mat)
  
  strain_coverage <- colSums(str_mat)

  return(list(species_richness = sp_richness, strain_coverage = strain_coverage))
}

# p1 = ggplot(data.frame(day = se_q@colData$days, alpha_per_sample = alpha$strain_coverage),
#             aes(x = day, y = alpha_per_sample, fill = day)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.7) +   # smoother look, hide outliers if needed
#   geom_jitter(width = 0.2, alpha = 0.5, size = 1) + # optional: show points
#   # facet_grid(~study) +
#   theme_minimal() +
#   labs(x = "day", y = "Strain Coverage") +
#   theme(legend.position = "none")
# 
# p2 = ggplot(data.frame(day = se_q@colData$days, alpha_per_sample = alpha$species_richness),
#             aes(x = day, y = alpha_per_sample, fill = day)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.7) +   # smoother look, hide outliers if needed
#   geom_jitter(width = 0.2, alpha = 0.5, size = 1) + # optional: show points
#   # facet_grid(~study) +
#   theme_minimal() +
#   labs(x = "day", y = "Species richness") +
#   theme(legend.position = "none")
# 
# 
# p1/p2
# 
# p3 = ggplot(data.frame(disease = sy@colData$tumour_stage_AJCC, study = sy@colData$study, alpha_per_sample = alpha$strain_coverage),
#             aes(x = disease, y = alpha_per_sample, fill = disease)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.7) +   # smoother look, hide outliers if needed
#   geom_jitter(width = 0.2, alpha = 0.5, size = 1) + # optional: show points
#   # facet_grid(~study) +
#   theme_minimal() +
#   labs(x = "Disease status", y = "Strain Coverage") +
#   theme(legend.position = "none")
# 
# p4 = ggplot(data.frame(disease = sy@colData$tumour_stage_AJCC, study = sy@colData$study, alpha_per_sample = alpha$species_richness),
#             aes(x = disease, y = alpha_per_sample, fill = disease)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.7) +   # smoother look, hide outliers if needed
#   geom_jitter(width = 0.2, alpha = 0.5, size = 1) + # optional: show points
#   # facet_grid(~study) +
#   theme_minimal() +
#   labs(x = "Disease status", y = "Species richness") +
#   theme(legend.position = "none")
# 
# 
# p3/p4

# 
# 
# p1 = ggplot(data.frame(disease = sy@colData$tumour_location, study = sy@colData$study, alpha_per_sample = alpha$strain_coverage),
#        aes(x = disease, y = alpha_per_sample, fill = disease)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.7) +   # smoother look, hide outliers if needed
#   geom_jitter(width = 0.2, alpha = 0.5, size = 1) + # optional: show points
#   # facet_grid(~study) +
#   theme_minimal() +
#   labs(x = "Disease status", y = "Strain Coverage") +
#   theme(legend.position = "none")
# 
# p2 = ggplot(data.frame(disease = sy@colData$tumour_location, study = sy@colData$study, alpha_per_sample = alpha$species_richness),
#        aes(x = disease, y = alpha_per_sample, fill = disease)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.7) +   # smoother look, hide outliers if needed
#   geom_jitter(width = 0.2, alpha = 0.5, size = 1) + # optional: show points
#   # facet_grid(~study) +
#   theme_minimal() +
#   labs(x = "Disease status", y = "Species richness") +
#   theme(legend.position = "none")
# 
# 
# p1/p2
# 
# p3 = ggplot(data.frame(disease = sy@colData$tumour_stage_AJCC, study = sy@colData$study, alpha_per_sample = alpha$strain_coverage),
#             aes(x = disease, y = alpha_per_sample, fill = disease)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.7) +   # smoother look, hide outliers if needed
#   geom_jitter(width = 0.2, alpha = 0.5, size = 1) + # optional: show points
#   # facet_grid(~study) +
#   theme_minimal() +
#   labs(x = "Disease status", y = "Strain Coverage") +
#   theme(legend.position = "none")
# 
# p4 = ggplot(data.frame(disease = sy@colData$tumour_stage_AJCC, study = sy@colData$study, alpha_per_sample = alpha$species_richness),
#             aes(x = disease, y = alpha_per_sample, fill = disease)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.7) +   # smoother look, hide outliers if needed
#   geom_jitter(width = 0.2, alpha = 0.5, size = 1) + # optional: show points
#   # facet_grid(~study) +
#   theme_minimal() +
#   labs(x = "Disease status", y = "Species richness") +
#   theme(legend.position = "none")
# 
# 
# p3/p4

# 
# 
# estimate_beta_diversity = function(){
#   
#   
#   
#   if (!requireNamespace("vegan", quietly = TRUE)) {
#     stop("The 'began' package is required but is not installed. Please install it first.")
#   }
#   
#   asy = t(as.matrix(SummarizedExperiment::assay(sy)))
#   asy[asy > 95] = 1
#   
#   richness <- vegan::specnumber(asy)
#   dist_jaccard <- vegan::vegdist((asy>0)+0, method = "jaccard")
#   
#   # asy_h <- decostand(asy, method = "hellinger")
#   # 
#   # jac_dist <- vegdist(asy, method = "jaccard", binary = TRUE)
#   # 
#   hc <- hclust(dist_jaccard, method = "average")
#   plot(hc)
#   
#   nmds <- metaMDS(dist_jaccard, k=2)
#   ord_df <- as.data.frame(scores(nmds))
#   ord_df$Group <-  SummarizedExperiment::colData(sy)$Case_status  # replace with your metadata column
#   
#   library(ggplot2)
#   ggplot(ord_df, aes(NMDS1, NMDS2, color=Group)) +
#     geom_point(size=3) +
#     stat_ellipse() +
#     theme_bw() +
#     labs(title="NMDS (Jaccard distance)")
#   
#   
#   adonis2(dist_jaccard ~ Group, data=ord_df)
#   
#   hc <- hclust(dist_jaccard)
#   plot(hc, labels=ord_df$Group)
#   
#   
# }
# 
# 
