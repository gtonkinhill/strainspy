#' @importFrom S4Vectors DataFrame

setClassUnion("DataFrameOrNULL", c("DataFrame", "NULL"))

#' strainspy_fit: A Class for storing fit data from all strainspy models
#'
#' The `strainspy_fit` class is designed to store the results of strainspy
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
#' @slot residuals A `DataFrame` containing residuals from the model.
#' @slot convergence A logical indicating whether model fitting converged.
#' @slot design A matrix representing the design matrix of the model.
#' @slot assay A matrix containing the assay data used for fitting the model.
#' @slot call The matched call of the model.
#'
#' @export
methods::setClass("strainspy_fit",
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

#' Show Method for strainspy_fit
#'
#' Displays a summary of a `strainspy_fit` object, including model design,
#' convergence, and coefficient matrix dimensions.
#'
#' @param object An object of class `strainspy_fit`.
#' @export
methods::setMethod("show", "strainspy_fit", function(object) {
  cat("strainspy_fit object\n")
  cat(rep("-", 30), "\n", sep = "")
  
  cat("Design formula: ")
  print(object@design)
  
  cat("Converged=", length(which(object@convergence==T)) , ", Failed=", length(which(object@convergence==F)), "\n", sep = '')
  
  cat("Coefficients: ", nrow(object@coefficients), " features x ", ncol(object@coefficients), " terms\n")
  
  if (!is.null(object@zi_coefficients)) {
    cat("ZI Coefficients: ", nrow(object@zi_coefficients), " features x ", ncol(object@zi_coefficients), " terms\n")
  } else {
    cat("ZI Coefficients: None\n")
  }
  
  cat("Residuals available: ", ifelse(!is.null(object@residuals), "Yes", "No"), "\n")
  cat("Call:\n")
  print(object@call)
})

#' Access Zero-Inflation Coefficients
#'
#' Returns the zero-inflation coefficients from a `strainspy_fit` object.
#'
#' @param object An object of class `strainspy_fit`.
#' @return A `DataFrame` containing zero-inflation coefficients.
#' @export
methods::setGeneric("getZICoefficients", function(object) methods::standardGeneric("getZICoefficients"))

#' @rdname getZICoefficients
#' @export
methods::setMethod("getZICoefficients", "strainspy_fit", function(object) {
  object@zi_coefficients
})

#' Access Residuals
#'
#' Returns the residuals from a `strainspy_fit` object.
#'
#' @param object An object of class `strainspy_fit`.
#' @return A `DataFrame` containing residuals.
#' @export
methods::setGeneric("getResiduals", function(object) methods::standardGeneric("getResiduals"))

#' @rdname getResiduals
#' @export
methods::setMethod("getResiduals", "strainspy_fit", function(object) {
  object@residuals
})

#' Access False Discovery Rates (FDR)
#'
#' Returns the false discovery rates (FDR) from a `strainspy_fit` object, if available.
#'
#' @param object An object of class `strainspy_fit`.
#' @return A `DataFrame` or `NULL` containing FDR values.
#' @export
methods::setGeneric("getFDR", function(object) methods::standardGeneric("getFDR"))

#' @rdname getFDR
#' @export
methods::setMethod("getFDR", "strainspy_fit", function(object) {
  object@fdr
})

#' Summary Method for strainspy_fit
#'
#' Provides a high-level summary of a `strainspy_fit` object, including its call,
#' coefficients, zero-inflation coefficients, dispersion, and the availability
#' of FDR and residuals.
#'
#' @param object An object of class `strainspy_fit`.
#' @return A list summarizing the `strainspy_fit` object.
#'
#' @importFrom utils head
#'
#' @export
methods::setMethod("summary", "strainspy_fit", function(object) {
  list(
    call = object@call,
    coefficients = head(object@coefficients),
    zi_coefficients = head(object@zi_coefficients),
    residuals_available = !is.null(object@residuals)
  )
})


#' Get Contig names in a fit object
#'
#' Gets the ordered list of fitted contig names in a `strainspy_fit`
#'
#' @param object An object of class `strainspy_fit`.
#' @return A vector of contig names
#' @export
methods::setGeneric("getContigNames", function(object) methods::standardGeneric("getContigNames"))

#' @rdname getContigNames
#' @export
methods::setMethod("getContigNames", "strainspy_fit", function(object) {
  object@row_data$Contig_name
})

#' Get Genomes in a fit object
#'
#' Gets the ordered list of fitted genomes in a `strainspy_fit`
#'
#' @param object An object of class `strainspy_fit`.
#' @return A vector of fitted genomes
#' @export
methods::setGeneric("getGenomes", function(object) methods::standardGeneric("getGenomes"))

#' @rdname getGenomes
#' @export
methods::setMethod("getGenomes", "strainspy_fit", function(object) {
  object@row_data$Genome_file
})

#' Get Sample names in a fit object
#'
#' Gets the ordered list of samples used in `strainspy_fit`
#'
#' @param object An object of class `strainspy_fit`.
#' @return A vector of fitted sample names 
#' @export
methods::setGeneric("getSamples", function(object) methods::standardGeneric("getSamples"))

#' @rdname getSamples
#' @export
methods::setMethod("getSamples", "strainspy_fit", function(object) {
  names(object@residuals)
})