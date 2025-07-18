#' strainspy_priors: A Class for storing model fitting priors
#'
#' The `strainspy_priors` class is designed to store MAP priors used for ordinal 
#' and zero-inflated beta model fitting using `glmmTMB`.
#'
#' @slot method Character. Prior estimation method used ("preset_weak", "preset_strong", "empirical", "manual").
#' @slot priors_df Data frame. Data frame containing prior parameters.
#' @slot boot_fixef Numeric matrix. Bootstrap SD values for fixed effect priors
#' @slot boot_fixef_ZI Numeric matrix. Bootstrap SD values for zero inflated fixed effect priors
#' @slot design The formula with with effect terms used
#' @slot call The matched call of the model.
#'
#' @export
methods::setClass("strainspy_priors",
                  slots = list(
                    method = "character",
                    priors_df = "data.frame",
                    boot_fixef = "matrix",
                    boot_fixef_ZI = "matrix",
                    design = "formula",
                    call = "call"
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
  invisible(x)
})

#' Show method for strainspy_priors objects
#'
#' @param object An object of class \code{strainspy_priors}
#' @return Invisibly returns the object
#' @export
methods::setMethod("show", "strainspy_priors", function(object) {
  cat("strainspy_priors object\n")
  cat("-------------------------\n")
  cat("Method: ", object@method, "\n\n")
  
  # Priors table summary
  if (!is.null(object@priors_df) && nrow(object@priors_df) > 0) {
    cat("Priors:\n")
    tab <- table(object@priors_df$class)
    for (cl in names(tab)) {
      cat(" -", cl, ":", tab[[cl]], "coefficients\n")
    }
    
    cat("\nExample priors:\n")
    print(utils::head(object@priors_df, 5))
    if (nrow(object@priors_df) > 5) {
      cat("... and", nrow(object@priors_df) - 5, "more rows\n")
    }
  } else {
    cat("No priors specified.\n")
  }
  
  # Bootstrap summaries
  warn_outlier_sd <- function(mat, label) {
    if (is.null(mat)) return()
    sds <- apply(mat, 1, stats::median)
    n_low <- sum(sds < 1)
    n_high <- sum(sds > 5)
    if (n_low > 0 || n_high > 0) {
      cat(sprintf("\nWarning: %d %s prior(s) may be too strong (SD < 1), and %d may be too weak (SD > 5).\n",
                  n_low, label, n_high))
    }
  }
  
  if (!is.null(object@boot_fixef)) {
    cat("\nBootstrap summary for fixef priors:\n")
    print(summary(as.vector(object@boot_fixef)))
    warn_outlier_sd(object@boot_fixef, "fixef")
  }
  
  if (!is.null(object@boot_fixef_ZI)) {
    cat("\nBootstrap summary for fixef_zi priors:\n")
    print(summary(as.vector(object@boot_fixef_ZI)))
    warn_outlier_sd(object@boot_fixef_ZI, "fixef_zi")
  }
  
  invisible(object)
})

#' Plot bootstrap prior SDs (fixef + fixef_zi) for a coefficient
#' 
#' @importFrom stats quantile
#'
#' @param object A \code{strainspy_priors} object.
#' @param prior Character. Coefficient name (must exist in both prior matrices).
#' @return A patchwork side-by-side ggplot object.
#' @rdname plot_prior_bootstrap
#' @export
methods::setGeneric("plot_prior_bootstrap", function(object, prior) standardGeneric("plot_prior_bootstrap"))

#' @rdname plot_prior_bootstrap
methods::setMethod("plot_prior_bootstrap", signature(object = "strainspy_priors"),
                   function(object, prior) {
                     
                     if(object@method != "empirical"){
                       stop("Bootstrap priors are only available for empirical bayes priors")
                     }
                     
                     if(length(object@boot_fixef) == 0 | length(object@boot_fixef_ZI) == 0){
                       stop("Bootstrap variables are empty")
                     }
                     
                     if(!all.equal(dim(object@boot_fixef_ZI), dim(object@boot_fixef))){
                       stop("Bootstrap variables have different dimensions")
                     }
                     
                     requireNamespace("ggplot2")
                     requireNamespace("patchwork")
                     
                     get_plot <- function(boot_sd, title) {
                       boot_sd <- as.numeric(boot_sd)
                       ci <- stats::quantile(boot_sd, c(0.025, 0.975))
                       mean_val <- mean(boot_sd)
                       median_val <- median(boot_sd)
                       
                       if (median_val < 1)
                         warning(sprintf("Estimated prior SD for %s (%s) is %.2f: this may be too strong and reduce sensitivity.",
                                         prior, title, median_val), call. = FALSE)
                       if (median_val > 5)
                         warning(sprintf("Estimated prior SD for %s (%s) is %.2f: this may be too weak and reduce specificity.",
                                         prior, title, median_val), call. = FALSE)
                       
                       df <- data.frame(sd = boot_sd)
                       
                       ggplot2::ggplot(df, ggplot2::aes(x = sd)) +
                         ggplot2::geom_histogram(bins = 50, fill = "grey80", color = "white") +
                         ggplot2::geom_vline(xintercept = mean_val, color = "blue", linetype = "dashed") +
                         ggplot2::geom_vline(xintercept = median_val, color = "red", linetype = "dashed") +
                         ggplot2::geom_vline(xintercept = ci, color = "darkgreen", linetype = "dotted") +
                         ggplot2::labs(
                           title = paste(title, "-", prior),
                           x = "Bootstrap SD values", y = "Frequency"
                         ) +
                         ggplot2::theme_minimal(base_size = 14)
                     }
                     
                     if (!(prior %in% rownames(object@boot_fixef)) || !(prior %in% rownames(object@boot_fixef_ZI))) {
                       stop("Prior not found in boot_fixef or boot_fixef_ZI")
                     }
                     
                     p1 <- get_plot(object@boot_fixef[prior, ], "fixef")
                     p2 <- get_plot(object@boot_fixef_ZI[prior, ], "fixef_ZI")
                     
                     p1 + p2 + patchwork::plot_layout(guides = "collect")
                   })


#' Access fixef priors needed for ordinal beta regression
#'
#' @param object A \code{strainspy_priors} object.
#' @return A data frame of prior parameters for fixed effects.
#' @export
methods::setGeneric("extract_fixef_priors", function(object) standardGeneric("extract_fixef_priors"))

#' @rdname extract_fixef_priors
#' @export
methods::setMethod("extract_fixef_priors", "strainspy_priors", function(object) {
  if (is.null(object)) {
    NULL
  } else {
    object@priors_df[object@priors_df$class == "fixef", ]
  }
})