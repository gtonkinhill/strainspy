#' strainspy_priors: A Class for storing model fitting priors
#'
#' The `strainspy_priors` class is designed to store MAP priors used for ordinal 
#' and zero-inflated beta model fitting using `GLMMTMB`.
#'
#' @slot method Character. Prior estimation method used ("preset_weak", "preset_strong", "empirical", "manual").
#' @slot priors_df Data frame. Data frame containing prior parameters.
#'
#' @export
methods::setClass("strainspy_priors",
                  slots = list(
                    method = "character",
                    priors_df = "data.frame"
                  )
)



#' Print method for strainspy_priors objects
#'
#' @param x An object of class \code{strainspy_priors}
#' @param ... Additional arguments (unused)
#' @return Prints a concise summary of the object
#' @export
methods::setMethod("print", "strainspy_priors", function(x, ...) {
  cat("strainspy_priors object\n")
  cat("Method used: ", x@method, "\n")
  if (nrow(x@priors_df) > 0) {
    cat("Number of prior entries: ", nrow(x@priors_df), "\n")
  } else {
    cat("No priors specified.\n")
  }
  if (nzchar(x@notes)) {
    cat("Notes:\n", x@notes, "\n")
  }
  invisible(x)
})

#' Access prior dataframe
#'
#' @param object strainspy_priors object
#' @return Data frame of prior parameters
#' @export
methods::setGeneric("getPriors", function(object) standardGeneric("getPriors"))

#' @rdname getPriors
#' @export
methods::setMethod("getPriors", "strainspy_priors", function(object) {
  object@priors_df
})


get_preset_priors = function(nbeta, type = c("weak", "strong")) {
  if (type == "weak") {
    # Example weak priors
    data.frame(
      prior = rep("normal(0,5)", 2*nbeta),
      class = rep(c("fixef", "fixef_zi"), each=nbeta),
      coef  = rep(as.character(seq(1,nbeta)), 2))
  } else {
    # Example strong priors
    data.frame(
      prior = rep("normal(0,1)", 2*nbeta),
      class = rep(c("fixef", "fixef_zi"), each=nbeta),
      coef  = rep(as.character(seq(1,nbeta)), 2))
  }
}