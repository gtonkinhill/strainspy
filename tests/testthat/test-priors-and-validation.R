testthat::test_that("define_priors supports preset and manual validation", {
  d <- ss_test_data()
  design <- stats::as.formula("~ Case_status + Age_at_collection")

  testthat::expect_message(
    strainspy::define_priors(d$se_filtered[1:20, ], design, method = "preset_weak"),
    "No MAP prior will be set on the intercept"
  )

  pri <- suppressMessages(
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

# --- class validity -----------------------------------------------------------
# A malformed priors_df draws only "replacement has length zero" from glmmTMB,
# and because the fitters swallow per-feature errors the whole fit comes back
# with zero rows. Validity turns that into one clear message up front.

mk_priors <- function(...) {
  methods::new("strainspy_priors", boot_fixef = matrix(0, 0, 0),
               boot_fixef_ZI = matrix(0, 0, 0), design = ~ x,
               call = call("f"), ...)
}

testthat::test_that("strainspy_priors rejects a priors_df glmmTMB cannot use", {
  testthat::expect_error(
    mk_priors(method = "manual", priors_df = data.frame(wrong_col = "normal(0,1)")),
    "missing the column"
  )
  # An empty priors_df is legitimate: no priors were set on any term.
  testthat::expect_s4_class(
    mk_priors(method = "manual", priors_df = data.frame()), "strainspy_priors"
  )
  testthat::expect_s4_class(
    mk_priors(method = "empirical",
              priors_df = data.frame(prior = "normal(0,1)", class = "fixef",
                                     coef = "2")),
    "strainspy_priors"
  )
})

testthat::test_that("strainspy_priors constrains the method slot", {
  testthat::expect_error(
    mk_priors(method = "nonsense", priors_df = data.frame()), "must be one of"
  )
  testthat::expect_error(
    mk_priors(method = c("a", "b"), priors_df = data.frame()),
    "single, non-missing"
  )
})

# --- resolve_priors -----------------------------------------------------------

testthat::test_that("resolve_priors treats only an explicit NULL as the MLE case", {
  # An explicit NULL means "fit without priors" and must keep working.
  testthat::expect_null(resolve_priors(NULL, NULL, NULL))

  # Anything unrecognised is a mistake, not a request for an unpenalised fit.
  # Returning NULL here would silently change the model that gets fitted.
  for (bad in list(42, TRUE, list(a = 1), factor("preset_weak"), matrix(1))) {
    testthat::expect_error(resolve_priors(bad, NULL, NULL), "must be NULL")
  }
})

testthat::test_that("resolve_priors validates character input", {
  testthat::expect_error(
    resolve_priors(c("preset_weak", "preset_strong"), NULL, NULL),
    "single character string"
  )
  testthat::expect_error(
    resolve_priors("not_a_preset", NULL, NULL), "must be one of"
  )
})
