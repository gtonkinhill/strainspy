testthat::test_that("top_hits and hadjust return expected tabular outputs", {
  d <- ss_test_data()
  fit_z <- ss_fit_zib_readme()

  th <- strainspy::top_hits(fit_z, coef = 2, alpha = 1)
  testthat::expect_s3_class(th, "tbl_df")
  testthat::expect_true(all(c("coefficient", "p_value", "p_adjust") %in% colnames(th)))
  testthat::expect_identical(attr(th, "phenotype_coef"), 2)

  hadj <- strainspy::hadjust(fit_z, coef = 2, taxonomy = d$taxonomy)
  testthat::expect_s3_class(hadj, "tbl_df")
  testthat::expect_true(all(c("Level", "Model", "Name", "p_adjust") %in% colnames(hadj)))
})

testthat::test_that("plot data-return modes are stable for downstream workflows", {
  d <- ss_test_data()
  fit_z <- ss_fit_zib_readme()

  volcano_df <- strainspy::plot_volcano(fit_z, coef = 2, plot = FALSE)
  testthat::expect_s3_class(volcano_df, "tbl_df")
  testthat::expect_true(all(c("Coefficient", "p", "p_adj", "Model") %in% colnames(volcano_df)))

  manhattan_df <- strainspy::plot_manhattan(fit_z, coef = 2, taxonomy = d$taxonomy, plot = FALSE)
  testthat::expect_s3_class(manhattan_df, "data.frame")
  testthat::expect_true(all(c("Level", "Model", "Name", "p_adjust", "index_min", "index_max") %in% colnames(manhattan_df)))

  testthat::expect_error(
    strainspy::plot_manhattan(fit_z, aggregate_by_taxa = TRUE, taxonomy = NULL),
    "Cannot aggregate by taxa"
  )
})

testthat::test_that("distribution plotting helpers return expected object types", {
  d <- ss_test_data()
  contigs <- rownames(d$se_filtered)[1:3]

  dist_df <- strainspy::plot_ani_dist(d$se_filtered, "Case_status", contigs, plot = FALSE)
  testthat::expect_s3_class(dist_df, "tbl_df")
  testthat::expect_true(all(c("sample", "Contig", "ANI", "Case_status") %in% colnames(dist_df)))

  hist_plot <- strainspy::plot_histogram(
    d$se_filtered,
    phenotype = "Case_status",
    contig = rownames(d$se_filtered)[1],
    bins = 40,
    fit_spline = FALSE,
    separate_facets = TRUE
  )
  testthat::expect_s3_class(hist_plot, "ggplot")

  summary_tbl <- attr(hist_plot, "summary_table")
  testthat::expect_s3_class(summary_tbl, "tbl_df")
  testthat::expect_true(all(c("stat", "value") %in% colnames(summary_tbl)))
})
