#' Provides the top results of a betaGLM object, sorted by adjusted p-value.
#'
#' This function performs multiple testing adjustment on the results of a betaGLM object
#' integrating taxonomic information. Various correction methods are available,
#' including Harmonic Mean P-value (HMP), Bonferroni and Benjamini-Hochberg (BH).
#'
#' @param object A `betaGLM` object.
#' @param coef The coefficient to report and use for p-value adjustment. Defaults to 2.
#' @param method The method for p-value adjustment. Available methods are a subset of those in `p.adjust`,
#'   specifically those valid for dependent tests. Common choices include "holm", "hochberg", "hommel",
#'   and "BH" (Benjamini-Hochberg)
#' @param alpha The significance level to use for filtering results. Defaults to 0.1.
#' @return A tibble with the top hits, sorted by adjusted p-value.
#' @export
top_hits <- function(object, coef=2, method = "holm", alpha=0.1) {
  # Validate input
  if (!inherits(object, "betaGLM")) {
    stop("Input must be a betaGLM object.")
  }

  # Check method is one of the available options
  if (!method %in% c("bonferroni", "BH", "BY", "holm")) {
    stop("Method must be one of 'bonferroni', 'BH', or 'holm'.")
  }

  # Check coef index validity
  if (coef > length(slot(object, "coefficients")) || coef < 1) {
    stop("Invalid coefficient index: ", coef)
  }

  # Extract p-values & coefficients
  res <- tibble::as_tibble(slot(object, 'row_data')) |>
    tibble::add_column(coefficient=slot(object, 'coefficients')[[coef]]) |>
    tibble::add_column(std_error=slot(object, 'std_errors')[[coef]]) |>
    tibble::add_column(p_value=slot(object, 'p_values')[[coef]]) |>
    tibble::add_column(p_adjust=p.adjust(slot(object, 'p_values')[[coef]], method = method)) |>
    dplyr::arrange(p_adjust)

  # Add zero-inflated coefficients if available
  if (!is.null(object@zi_coefficients)){
    res <- res |>
      tibble::add_column(zi_coefficient=object@zi_coefficients[[coef]]) |>
      tibble::add_column(zi_std_error=object@zi_std_errors[[coef]]) |>
      tibble::add_column(zi_p_value=object@zi_p_values[[coef]]) |>
      tibble::add_column(zi_p_adjust=p.adjust(object@zi_p_values[[coef]], method = method)) |>
      dplyr::arrange(pmin(p_adjust, zi_p_adjust))
  }

  # Filter on provided alpha
  res <- res |> dplyr::filter(pmin(p_adjust, zi_p_adjust) <= alpha)

  return(res)
}
