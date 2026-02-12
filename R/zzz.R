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
  "zero_count", "total", "zero_prop", "Value_orig", "value"
))


## Quick sizeable data test set for use: - KEEP COMMENTED
# meta_path <- "../strainspy-manuscript/data/ash_pancancer/metadata_full.tsv"
# meta <- read.csv(meta_path, sep = '\t')
# meta = cbind(run_acc = meta$run_accession, meta)
# meta$RvsP = "R"
# meta$RvsP[which(meta$BOR == "PD" | meta$BOR == "cPD")] = "NR"
# meta$RvsP = factor(meta$RvsP, levels = c("NR", "R"))
# 
# sy = read_sylph("../strainspy-manuscript/data/ash_pancancer/combined_q_99.tsv.gz")
# colnames(sy) <- gsub("_1", "", colnames(sy))
# SummarizedExperiment::colData(sy)$Sample_file <- gsub("_1", "", basename(SummarizedExperiment::colData(sy)$Sample_file))
# 
# meta = meta[match(colnames(sy), meta$run_accession), ]
# rmidx = which(meta$BOR == "SD")
# if(length(rmidx) > 0){
#   meta = meta[-rmidx, ]
#   sy = sy[, -rmidx]
# }
# 
# sy <- filter_by_presence(sy, min_nonzero = 8)
# all(colnames(sy) %in% meta$run_accession)
# all(meta$run_accession %in% colnames(sy))
# 
# sy = modify_metadata(sy, meta)
# design <- as.formula("~ RvsP + histology_cohort.x") # age + sex + BMI + ECOG_baseline + chemo_use + antibiotic_use + PPI_use +  + NLR + PLAT + ALB + LDH")
