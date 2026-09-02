# The fitting backends signal a failed per-feature fit in two different ways, and
# `failed_fit_index()` has to recognise both without mistaking a backend that
# simply does not report a log-likelihood for one whose every fit failed.

ok_ll <- function(v) list(coefficients = matrix(1, 2, 4), residuals = 1:3,
                          log_likelihood = v)
no_ll <- function() list(coefficients = matrix(1, 2, 4), residuals = 1:3)

testthat::test_that("failed_fit_index() detects the NA log-likelihood protocol", {
  # glmZiBFit returns a fully-populated, NA-filled result for a failed fit, so
  # the zero-length test never fires and the NA log-likelihood is the only signal.
  res <- list(ok_ll(-5), ok_ll(NA), ok_ll(-3))
  testthat::expect_identical(failed_fit_index(res), 2L)
  testthat::expect_identical(failed_fit_index(list(ok_ll(NaN), ok_ll(-3))), 1L)
  testthat::expect_identical(failed_fit_index(list(ok_ll(NA), ok_ll(NA))), c(1L, 2L))
})

testthat::test_that("failed_fit_index() detects the NULL protocol", {
  # glmObFit, glmQBFit and abundanceFit return NULL for a failed fit.
  testthat::expect_identical(failed_fit_index(list(ok_ll(-5), NULL, ok_ll(-3))), 2L)
  testthat::expect_identical(failed_fit_index(list(no_ll(), no_ll(), NULL)), 3L)
})

testthat::test_that("a backend that never reports a log-likelihood loses no rows", {
  # Regression guard: glmQBFit has log_likelihood commented out of its return
  # list, so testing it blindly would flag every successful fit as failed.
  testthat::expect_identical(failed_fit_index(list(no_ll(), no_ll(), no_ll())),
                             integer(0))
})

testthat::test_that("failed_fit_index() stays aligned with the results list", {
  # unlist() drops the NULL log-likelihood of a zero-length entry, shortening the
  # vector and shifting every later index. Positions must refer to `results`.
  res <- list(NULL, ok_ll(NA), ok_ll(-5), ok_ll(-3))
  testthat::expect_identical(failed_fit_index(res), c(1L, 2L))
  testthat::expect_length(unlist(lapply(res, function(x) x$log_likelihood)), 3L)
})

testthat::test_that("failed_fit_index() handles edge cases", {
  testthat::expect_identical(failed_fit_index(list()), integer(0))
  testthat::expect_identical(failed_fit_index(list(ok_ll(-5), ok_ll(-3))), integer(0))
  # A malformed, non-scalar log-likelihood counts as a failure.
  testthat::expect_identical(failed_fit_index(list(ok_ll(c(-1, -2)), ok_ll(-3))), 1L)
})
