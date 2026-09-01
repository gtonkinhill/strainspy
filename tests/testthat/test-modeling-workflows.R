testthat::test_that("glmZiBFit README-like workflow returns valid fit object", {
  fit_z <- ss_fit_zib()
  fit_readme <- ss_fit_zib_readme()

  testthat::expect_s4_class(fit_z, "strainspy_fit")
  testthat::expect_s4_class(fit_readme, "strainspy_fit")

  testthat::expect_equal(nrow(fit_z@coefficients), nrow(fit_z@row_data))
  testthat::expect_equal(nrow(fit_readme@coefficients), nrow(fit_readme@row_data))

  testthat::expect_false(is.null(fit_z@zi_coefficients))
  testthat::expect_false(is.null(fit_readme@zi_coefficients))

  testthat::expect_true(is.logical(fit_z@convergence))
  testthat::expect_true(is.logical(fit_readme@convergence))
})

testthat::test_that("glmQBFit remains available as secondary model", {
  fit_q <- ss_fit_qb()

  testthat::expect_s4_class(fit_q, "strainspy_fit")
  testthat::expect_equal(nrow(fit_q@coefficients), nrow(fit_q@row_data))
  testthat::expect_null(fit_q@zi_coefficients)
})

testthat::test_that("glmObFit and abundanceFit run on small examples", {
  d <- ss_test_data()
  design <- stats::as.formula("~ Case_status + Age_at_collection")

  fit_ob <- suppressWarnings(
    strainspy::glmObFit(d$se_filtered[1:15, ], design, nthreads = 1)
  )
  testthat::expect_s4_class(fit_ob, "strainspy_fit")
  testthat::expect_equal(nrow(fit_ob@coefficients), nrow(fit_ob@row_data))

  se_abund <- strainspy::read_sylph(
    d$sylph_path,
    meta_data = d$meta,
    variable = "Taxonomic_abundance"
  )
  fit_abund <- strainspy::abundanceFit(se_abund[1:20, ], design, nthreads = 1)

  testthat::expect_s4_class(fit_abund, "strainspy_fit")
  testthat::expect_null(fit_abund@zi_coefficients)
  testthat::expect_equal(nrow(fit_abund@coefficients), nrow(fit_abund@row_data))
})

testthat::test_that("caseControlFit runs and enforces formula requirements", {
  d <- ss_test_data()

  design_ok <- stats::as.formula("Case_status ~ Value + Age_at_collection")
  fit_cc <- strainspy::caseControlFit(d$se_filtered[1:12, ], design_ok, nthreads = 1)

  testthat::expect_s4_class(fit_cc, "strainspy_fit")
  testthat::expect_equal(nrow(fit_cc@coefficients), nrow(fit_cc@row_data))

  design_bad <- stats::as.formula("Case_status ~ Age_at_collection")
  testthat::expect_error(
    strainspy::caseControlFit(d$se_filtered[1:12, ], design_bad, nthreads = 1),
    "must contain the term Value"
  )
})
