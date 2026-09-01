testthat::test_that("define_priors supports preset and manual validation", {
  d <- ss_test_data()
  design <- stats::as.formula("~ Case_status + Age_at_collection")

  testthat::expect_warning(
    strainspy::define_priors(d$se_filtered[1:20, ], design, method = "preset_weak"),
    "No MAP prior will be set to intercept"
  )

  pri <- suppressWarnings(
    strainspy::define_priors(d$se_filtered[1:20, ], design, method = "preset_weak")
  )

  testthat::expect_s4_class(pri, "strainspy_priors")
  testthat::expect_true(nrow(pri@priors_df) > 0)

  manual_df <- data.frame(
    prior = c("normal(0,2)", "normal(0,2)"),
    class = c("fixef", "fixef_zi"),
    coef = c("2", "2"),
    stringsAsFactors = FALSE
  )

  testthat::expect_warning(
    strainspy::define_priors(
      d$se_filtered[1:20, ],
      design,
      method = "manual",
      priors_df = manual_df
    ),
    "No checks will be performed"
  )

  pri_manual <- suppressWarnings(
    strainspy::define_priors(
      d$se_filtered[1:20, ],
      design,
      method = "manual",
      priors_df = manual_df
    )
  )

  testthat::expect_s4_class(pri_manual, "strainspy_priors")

  testthat::expect_error(
    strainspy::define_priors(d$se_filtered[1:20, ], design, method = "manual", priors_df = data.frame(prior = "normal(0,1)")),
    "missing required column"
  )

  testthat::expect_error(
    strainspy::define_priors(d$se_filtered[1:20, ], design, method = "unknown_method"),
    "Unknown method"
  )
})

testthat::test_that("top_hits validates coef, alpha, and method inputs", {
  fit <- ss_fit_zib()

  testthat::expect_error(strainspy::top_hits(fit, coef = 0), "out of bounds")
  testthat::expect_error(strainspy::top_hits(fit, coef = 2.5), "single integer index")
  testthat::expect_error(strainspy::top_hits(fit, alpha = 1.1), "between 0 and 1")
  testthat::expect_error(strainspy::top_hits(fit, method = "hochberg"), "Method must be one of")
})

testthat::test_that("glmZiBFit validates method argument", {
  d <- ss_test_data()
  design <- stats::as.formula("~ Case_status + Age_at_collection")

  testthat::expect_error(
    strainspy::glmZiBFit(d$se_filtered[1:10, ], design, method = "not_a_method", nthreads = 1),
    "should be one of"
  )
})
