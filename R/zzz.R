# Non-standard-evaluation column and aesthetic names referenced inside
# data.table, dplyr and ggplot2 expressions. Function names do NOT belong
# here - declare those with @importFrom in R/strainspy-package.R instead, so
# that R CMD check can still catch a genuinely missing import.
utils::globalVariables(c(
  ".", "#clade_name", "ANI", "Class",
  "coefficient", "Coefficient", "col_indices", "comp_nz",
  "comp_total", "component", "component_p", "Contig",
  "Contig_name", "density", "Dim1", "Dim2",
  "DistanceToCentroid", "Domain", "estimate", "fail_reasons",
  "Family", "genome", "Genome", "Genome_file",
  "Genus", "Group", "index", "index_max",
  "index_min", "isTip", "label", "Level",
  "log_p_adjust", "log10p", "max_index", "min_index",
  "Model", "Name", "NMDS1", "NMDS2",
  "Order", "p_adjust", "p_value", "PC1",
  "PC2", "phenotype", "Phylum", "pred",
  "query_name", "ref_nz", "ref_total", "relative_abundance",
  "row_indices", "Sample_file", "Species", "spline_density",
  "spline_x", "tax", "taxonomy", "Taxonomy",
  "total", "V1", "V2", "value",
  "Value_orig", "y", "zero_count", "zero_prop",
  "zi_coefficient", "zi_p_adjust", "zi_p_value"
))


## Quick sizeable data test set for use: - KEEP COMMENTED
# meta_path <- "../strainspy-manuscript/data/ash_pancancer/metadata_full.tsv"
# meta <- read.csv(meta_path, sep = '\t')
# meta = cbind(run_acc = meta$run_accession, meta)
# meta$RvsP = "R"
# meta$RvsP[which(meta$BOR == "PD" | meta$BOR == "cPD")] = "NR"
# meta$RvsP = factor(meta$RvsP, levels = c("NR", "R"))
# 
# se = read_sylph("../strainspy-manuscript/data/ash_pancancer/combined_q_99.tsv.gz")
# colnames(se) <- gsub("_1", "", colnames(se))
# SummarizedExperiment::colData(se)$Sample_file <- gsub("_1", "", basename(SummarizedExperiment::colData(se)$Sample_file))
# 
# meta = meta[match(colnames(se), meta$run_accession), ]
# rmidx = which(meta$BOR == "SD")
# if(length(rmidx) > 0){
#   meta = meta[-rmidx, ]
#   se = se[, -rmidx]
# }
# 
# se <- filter_by_presence(se, min_nonzero = 8)
# all(colnames(se) %in% meta$run_accession)
# all(meta$run_accession %in% colnames(se))
# 
# se = modify_metadata(se, meta)
# design <- as.formula("~ RvsP + histology_cohort.x + age + sex + BMI") # age + sex + BMI + ECOG_baseline + chemo_use + antibiotic_use + PPI_use +  + NLR + PLAT + ALB + LDH")
