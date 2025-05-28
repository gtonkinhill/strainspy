#' @importFrom S4Vectors DataFrame

setClassUnion("DataFrameOrNULL", c("DataFrame", "NULL"))

#' betaGLM: A Class for Generalized Linear Models with Zero-Inflation
#'
#' The `betaGLM` class is designed to store the results of strainspy
#' 
#' @importClassesFrom S4Vectors DataFrame
#' @importMethodsFrom S4Vectors show
#'
#'
#' @slot row_data A `DataFrame` containing row data.
#' @slot coefficients A `DataFrame` containing model coefficients.
#' @slot std_errors A `DataFrame` containing standard errors of the coefficients.
#' @slot p_values A `DataFrame` containing p-values for the coefficients.
#' @slot zi_coefficients A `DataFrame` containing zero-inflation coefficients.
#' @slot zi_std_errors A `DataFrame` containing zero-inflation standard errors.
#' @slot zi_p_values A `DataFrame` containing zero-inflation p-values.
#' @slot fdr A `DataFrame` or `NULL`, containing false discovery rates (FDR).
#' @slot zi_fdr A `DataFrame` or `NULL`, containing zero-inflation FDR values.
#' @slot dispersion A numeric vector containing dispersion estimates for the model.
#' @slot residuals A `DataFrame` containing residuals from the model.
#' @slot design A matrix representing the design matrix of the model.
#' @slot assay A matrix containing the assay data used for fitting the model.
#' @slot call The matched call of the model.
#'
#' @export
methods::setClass("betaGLM",
                  slots = list(
                    row_data = "DataFrame",
                    coefficients = "DataFrame",
                    std_errors = "DataFrame",
                    p_values = "DataFrame",
                    zi_coefficients = "DataFrameOrNULL",
                    zi_std_errors = "DataFrameOrNULL",
                    zi_p_values = "DataFrameOrNULL",
                    residuals = "DataFrameOrNULL",
                    convergence = "logical",
                    design = "formula",
                    assay = "matrix",
                    call = "call"
                  ))

#' Show Method for betaGLM
#'
#' Displays a summary of a `betaGLM` object, including the number of coefficients.
#'
#' @param object An object of class `betaGLM`.
#' @export
methods::setMethod("show", "betaGLM", function(object) {
  cat("betaGLM object\n")
  cat("Number of coefficients:", nrow(object@coefficients), "\n")
  cat("Number of zero-inflation coefficients:", nrow(object@zi_coefficients), "\n")
  cat("Residuals available:", !is.null(object@residuals), "\n")
})

#' Access Zero-Inflation Coefficients
#'
#' Returns the zero-inflation coefficients from a `betaGLM` object.
#'
#' @param object An object of class `betaGLM`.
#' @return A `DataFrame` containing zero-inflation coefficients.
#' @export
methods::setGeneric("getZICoefficients", function(object) methods::standardGeneric("getZICoefficients"))
methods::setMethod("getZICoefficients", "betaGLM", function(object) {
  object@zi_coefficients
})

#' Access Residuals
#'
#' Returns the residuals from a `betaGLM` object.
#'
#' @param object An object of class `betaGLM`.
#' @return A `DataFrame` containing residuals.
#' @export
methods::setGeneric("getResiduals", function(object) methods::standardGeneric("getResiduals"))
methods::setMethod("getResiduals", "betaGLM", function(object) {
  object@residuals
})

#' Access False Discovery Rates (FDR)
#'
#' Returns the false discovery rates (FDR) from a `betaGLM` object, if available.
#'
#' @param object An object of class `betaGLM`.
#' @return A `DataFrame` or `NULL` containing FDR values.
#' @export
methods::setGeneric("getFDR", function(object) methods::standardGeneric("getFDR"))
methods::setMethod("getFDR", "betaGLM", function(object) {
  object@fdr
})

#' Summary Method for betaGLM
#'
#' Provides a high-level summary of a `betaGLM` object, including its call,
#' coefficients, zero-inflation coefficients, dispersion, and the availability
#' of FDR and residuals.
#'
#' @param object An object of class `betaGLM`.
#' @return A list summarizing the `betaGLM` object.
#'
#' @importFrom utils head
#'
#' @export
methods::setMethod("summary", "betaGLM", function(object) {
  list(
    call = object@call,
    coefficients = head(object@coefficients),
    zi_coefficients = head(object@zi_coefficients),
    residuals_available = !is.null(object@residuals)
  )
})
