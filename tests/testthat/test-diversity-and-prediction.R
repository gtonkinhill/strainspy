testthat::test_that("prep_for_prediction returns a model-ready data frame", {
  d <- ss_test_data()
  contigs <- rownames(d$se_filtered)[1:3]

  pred <- strainspy::prep_for_prediction(
    sy = d$se_filtered,
    outcome = "Case_status",
    contigs = contigs,
    covariates = "Age_at_collection",
    use_genome_names = FALSE
  )

  testthat::expect_s3_class(pred, "data.frame")
  testthat::expect_true("Case_status" %in% colnames(pred))
  testthat::expect_true(all(contigs %in% colnames(pred)))
  testthat::expect_equal(nrow(pred), ncol(d$se_filtered))

  testthat::expect_error(
    strainspy::prep_for_prediction(d$se_filtered, outcome = "MissingOutcome", contigs = contigs),
    "Outcome variable"
  )
})

testthat::test_that("richness and beta-diversity utilities return expected structures", {
  d <- ss_test_data()

  rich <- strainspy::estimate_sample_richness(d$se_filtered, d$taxonomy)
  testthat::expect_type(rich, "list")
  testthat::expect_true(all(c("species_richness", "strain_coverage") %in% names(rich)))
  testthat::expect_equal(length(rich$species_richness), ncol(d$se_filtered))

  dist_obj <- strainspy::estimate_beta_diversity(d$se_filtered, distance_only = TRUE)
  testthat::expect_s3_class(dist_obj, "dist")

  if (requireNamespace("vegan", quietly = TRUE)) {
    div_out <- suppressWarnings(
      strainspy::estimate_beta_diversity(
        d$se_filtered,
        phenotype = "Case_status",
        distance_only = FALSE,
        return_plots = TRUE
      )
    )
    testthat::expect_type(div_out, "list")
    testthat::expect_true(all(c("distance", "nmds", "ordination", "permanova", "betadisper", "plots") %in% names(div_out)))
    testthat::expect_true(all(c("NMDS", "BetaDispersionBox") %in% names(div_out$plots)))
  }
})

testthat::test_that("strainspy_fit accessor methods return coherent values", {
  fit_z <- ss_fit_zib()

  testthat::expect_equal(length(strainspy::getContigNames(fit_z)), nrow(fit_z@row_data))
  testthat::expect_equal(length(strainspy::getGenomes(fit_z)), nrow(fit_z@row_data))
  testthat::expect_s4_class(strainspy::getResiduals(fit_z), "DFrame")
  testthat::expect_s4_class(strainspy::getZICoefficients(fit_z), "DFrame")
})
