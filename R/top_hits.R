#' Provides the top results of a strainspy_fit object, sorted by adjusted p-value.
#'
#' This function performs multiple testing adjustment on the results of a strainspy_fit object
#' integrating taxonomic information. Various correction methods are available,
#' including Harmonic Mean P-value (HMP), Bonferroni and Benjamini-Hochberg (BH).
#'
#' @param object A `strainspy_fit` object.
#' @param coef The coefficient to report and use for p-value adjustment. Defaults to 2.
#' @param method The method for p-value adjustment. Available methods are a subset of those in `p.adjust`,
#'   specifically those valid for dependent tests. Common choices include "holm", "hochberg", "hommel",
#'   and "BH" (Benjamini-Hochberg) and "BY". 
#' @param alpha Numeric. Significance threshold for adjusted p-values. Defaults to 0.05.
#' 
#' @return A tibble with the top hits, sorted by adjusted p-value.
#' @export
#' @importFrom dplyr mutate
#' @importFrom tibble as_tibble add_column
top_hits <- function(object, coef=2, method = "holm", alpha=0.05) {
  # Validate input
  if (!inherits(object, "strainspy_fit")) {
    stop("Input must be a strainspy_fit")
  }

  # Check method is one of the available options
  if (!method %in% c("bonferroni", "BH", "BY", "holm")) {
    stop("Method must be one of 'bonferroni', 'BH', 'holm', 'BY', or 'holm'.")
  }

  # check if taxonomy is a string
  # if (!is.null(taxonomy) && is.character(taxonomy)) {
  #   tax <- read_taxonomy(taxonomy)
  # }

  # Extract p-values & coefficients
  res <- tibble::as_tibble(object@row_data) |>
    tibble::add_column(coefficient=object@coefficients[[coef]]) |>
    tibble::add_column(std_error=object@std_errors[[coef]]) |>
    tibble::add_column(p_value=object@p_values[[coef]]) |>
    tibble::add_column(p_adjust=p.adjust(object@p_values[[coef]], method = method)) # |>
    # dplyr::arrange(p_adjust)

  # Add zero-inflated coefficients if available
  if (!is.null(object@zi_coefficients)){
    res <- res |>
      tibble::add_column(zi_coefficient=object@zi_coefficients[[coef]]) |>
      tibble::add_column(zi_std_error=object@zi_std_errors[[coef]]) |>
      tibble::add_column(zi_p_value=object@zi_p_values[[coef]]) |>
      tibble::add_column(zi_p_adjust=p.adjust(object@zi_p_values[[coef]], method = method))  #|>
      # dplyr::arrange(pmin(p_adjust, zi_p_adjust))
  }
  
  ## Arrange by adjusted p_value
  if(!is.null(object@zi_coefficients)) {
    res <- res |> dplyr::arrange(pmin(p_adjust, zi_p_adjust))
  } else {
    res <- res |> dplyr::arrange(p_adjust)
  }

  # Filter on provided alpha
  if (!is.null(object@zi_coefficients)){
    res <- res |> dplyr::filter(pmin(p_adjust, zi_p_adjust) <= alpha)
  } else {
    res <- res |> dplyr::filter(p_adjust <= alpha)
  }

  return(res)
}
