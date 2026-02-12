#' Provides the top results of a strainspy_fit object, sorted by adjusted p-value.
#'
#' This function performs multiple testing adjustment on the results of a strainspy_fit object
#' integrating taxonomic information. Various correction methods are available,
#' including Harmonic Mean P-value (HMP), Bonferroni and Benjamini-Hochberg (BH).
#'
#' @param fit A `strainspy_fit` object.
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
top_hits <- function(fit, coef=2, method = "holm", alpha=0.05) {
  # Validate input
  if (!inherits(fit, "strainspy_fit")) {
    stop("Input must be a strainspy_fit")
  }
  
  # Check method is one of the available options
  if (!method %in% c("bonferroni", "BH", "BY", "holm")) {
    stop("Method must be one of 'bonferroni', 'BH', 'BY', or 'holm'.")
  }
  
 
  
  # Extract p-values & coefficients
  res <- tibble::as_tibble(fit@row_data) |>
    tibble::add_column(coefficient=fit@coefficients[[coef]]) |>
    tibble::add_column(std_error=fit@std_errors[[coef]]) |>
    tibble::add_column(p_value=fit@p_values[[coef]]) |>
    tibble::add_column(p_adjust=p.adjust(fit@p_values[[coef]], method = method)) # |>
  # dplyr::arrange(p_adjust)
  
  # Add zero-inflated coefficients if available
  if (!is.null(fit@zi_coefficients)){
    res <- res |>
      tibble::add_column(zi_coefficient=fit@zi_coefficients[[coef]]) |>
      tibble::add_column(zi_std_error=fit@zi_std_errors[[coef]]) |>
      tibble::add_column(zi_p_value=fit@zi_p_values[[coef]]) |>
      tibble::add_column(zi_p_adjust=p.adjust(fit@zi_p_values[[coef]], method = method))  #|>
    # dplyr::arrange(pmin(p_adjust, zi_p_adjust))
  }
  
  ## Arrange by adjusted p_value
  if(!is.null(fit@zi_coefficients)) {
    res <- res |> dplyr::arrange(pmin(p_adjust, zi_p_adjust))
  } else {
    res <- res |> dplyr::arrange(p_adjust)
  }
  
  # Filter on provided alpha
  if (!is.null(fit@zi_coefficients)){
    res <- res |> dplyr::filter(pmin(p_adjust, zi_p_adjust) <= alpha)
  } else {
    res <- res |> dplyr::filter(p_adjust <= alpha)
  }
  
  if(nrow(res) == 0) {
    warning(sprintf("Multiple testing correction using `%s`: No significant associations detected for coef = %d at alpha = %f", method, coef, alpha))
  } else {
    cat(paste("Found", nrow(res), "tophits for", names(fit@coefficients)[coef], "at alpha =", alpha, "using", method, "\n"))
    # shortcut to attach phenotype and other info to top_hits
    attr(res, "phenotype_coef") = coef
    # attr(res, "phenotype") = names(fit@coefficients)[coef] # This fails with multifactor variables
    attr(res, "method") = method
    attr(res, "alpha") = alpha
  }
  
  
  
  return(res)
}
