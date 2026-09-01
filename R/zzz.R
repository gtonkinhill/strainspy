utils::globalVariables(c(
  ".", "#clade_name", "ANI", "Class", "Coefficient", "Contig",
  "Contig_name", "Dim1", "Dim2", "Domain", "Family", "Genome_file", "Genus",
  "Level", "Model", "Name", "Order", "PC1", "PC2", "Phylum", "Sample_file",
  "Species", "Taxonomy", "V1", "V2", "abline", "ani", "as", "binomial",
  "capture.output", "coefficient", "col_indices", "dist", "genome", "index",
  "index_max", "index_min", "log10p", "log_p_adjust", "max_index", "min_index",
  "p.adjust", "p_adjust", "p_value", "phenotype", "prcomp", "pred", "query_name",
  "rbinom", "relative_abundance", "row_indices", "set", "slotNames",
  "str_extract", "strip_random_effects", "tax", "taxonomy", "tibble",
  "zi_coefficient", "zi_p_adjust", "zi_p_value", 'isTip', 'y', 'label', 'Genome',
  "zero_count", "total", "zero_prop", "Value_orig", "value", "component", 
  "ref_nz", "comp_nz", "ref_total", "comp_total", "estimate", "fail_reasons",
  "component_p", "NMDS1", "NMDS2", "Group", "DistanceToCentroid",
  "spline_x", "spline_density", "spline_y", "density"
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
