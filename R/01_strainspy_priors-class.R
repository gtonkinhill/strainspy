#' strainspy_priors: A Class for storing model fitting priors
#'
#' The `strainspy_priors` class is designed to store MAP priors used for ordinal 
#' and zero-inflated beta model fitting using `glmmTMB`.
#'
#' @slot method Character. Prior estimation method used ("preset_weak", "preset_strong", "empirical", "manual").
#' @slot priors_df Data frame. Data frame containing prior parameters.
#' @slot boot_fixef Numeric matrix. Bootstrap SD values for fixed effect priors
#' @slot boot_fixef_ZI Numeric matrix. Bootstrap SD values for zero inflated fixed effect priors
#' @slot boot_disp Numeric matrix. Bootstrap mean and SD values for the dispersion term of beta component.
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
                    boot_disp="matrix",
                    design = "formula",
                    call = "call"
                  )
)

# A malformed `priors_df` is not rejected by glmmTMB with anything useful: it
# fails with "replacement has length zero" deep inside the fit, and because the
# fitters wrap each feature in tryCatch(error = function(e) NULL), every feature
# is silently dropped and the user gets a fit object with zero rows. Validate
# the shape here so the problem is reported once, clearly, up front.
methods::setValidity("strainspy_priors", function(object) {
  msg <- character(0)

  valid_methods <- c("preset_weak", "preset_strong", "empirical", "manual")
  if (length(object@method) != 1L || is.na(object@method)) {
    msg <- c(msg, "`method` must be a single, non-missing character string.")
  } else if (!object@method %in% valid_methods) {
    msg <- c(msg, paste0("`method` must be one of ",
                         paste(sQuote(valid_methods), collapse = ", "),
                         "; got ", sQuote(object@method), "."))
  }

  # glmmTMB's `priors` argument requires these columns. An empty data frame is
  # legitimate: it means no priors were set on any term.
  if (nrow(object@priors_df) > 0L) {
    required <- c("prior", "class", "coef")
    missing_cols <- setdiff(required, colnames(object@priors_df))
    if (length(missing_cols) > 0L) {
      msg <- c(msg, paste0(
        "`priors_df` is missing the column(s) ",
        paste(sQuote(missing_cols), collapse = ", "),
        " required by `glmmTMB::glmmTMB(priors = )`. Present: ",
        paste(sQuote(colnames(object@priors_df)), collapse = ", "), "."
      ))
    }
  }

  if (length(msg) > 0L) msg else TRUE
})



#' Show method for strainspy_priors objects
#'
#' @param object An object of class \code{strainspy_priors}
#' @return Invisibly returns the object
#' @export
methods::setMethod("show", "strainspy_priors", function(object) {
  cat("strainspy_priors object\n")
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
  }
  if (!is.null(object@boot_fixef_ZI)) {
    cat("\nBootstrap summary for fixef_zi priors:\n")
    print(summary(as.vector(object@boot_fixef_ZI)))
  }
  
  warn_outlier_sd(object@boot_fixef, "fixef")
  warn_outlier_sd(object@boot_fixef_ZI, "fixef_zi")
  
  invisible(object)
})

#' Plot method for strainspy_priors objects
#'
#' @param object An object of class \code{strainspy_priors} that has bootstrap information
#' @param prior Prior to be plotted
#' @return Invisibly returns the object
#' @examples
#' methods::hasMethod("plot_prior_bootstrap", "strainspy_priors")
#' @export
methods::setGeneric("plot_prior_bootstrap", function(object, prior) standardGeneric("plot_prior_bootstrap"))

#' @rdname plot_prior_bootstrap
#' @export
methods::setMethod("plot_prior_bootstrap", signature(object = "strainspy_priors"),
                   function(object, prior) {
                     
                     if (object@method != "empirical") {
                       stop("Bootstrap plots are only available for 'empirical' method objects.")
                     }
                     
                     if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("patchwork", quietly = TRUE)) {
                       stop("Please install 'ggplot2' and 'patchwork' to use this plotting function.")
                     }
                     
                     # Internal plotting logic to keep code DRY
                     internal_plot <- function(vec, title) {
                       df <- data.frame(sd = as.numeric(vec))
                       med_val <- stats::median(df$sd)
                       
                       ggplot2::ggplot(df, ggplot2::aes(x = sd)) +
                         ggplot2::geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.8) +
                         ggplot2::geom_vline(xintercept = med_val, color = "red", linetype = "dashed") +
                         ggplot2::labs(title = title, subtitle = paste("Median SD:", round(med_val, 3)),
                                       x = "Bootstrap SD", y = "Count") +
                         ggplot2::theme_minimal()
                     }
                     
                     plots <- list()
                     
                     # Check fixef
                     if (prior %in% rownames(object@boot_fixef)) {
                       plots$p1 <- internal_plot(object@boot_fixef[prior, ], paste("Fixef:", prior))
                     }
                     
                     # Check ZI
                     if (prior %in% rownames(object@boot_fixef_ZI)) {
                       plots$p2 <- internal_plot(object@boot_fixef_ZI[prior, ], paste("Fixef_ZI:", prior))
                     }
                     
                     # Check Dispersion
                     if (prior == "dispersion" && length(object@boot_disp) > 0) {
                       plots$p1 <- internal_plot(object@boot_disp[1, ], "Dispersion (Term 1)")
                       if (nrow(object@boot_disp) > 1) {
                         plots$p2 <- internal_plot(object@boot_disp[2, ], "Dispersion (Term 2)")
                       }
                     }
                     
                     if (length(plots) == 0) stop("Coefficient ", prior, " not found in bootstrap matrices.")
                     
                     # Combine using patchwork
                     combined <- patchwork::wrap_plots(plots) + patchwork::plot_layout(guides = "collect")
                     return(combined)
                   }
)
#' Access fixef priors needed for ordinal beta regression
#'
#' @param object A \code{strainspy_priors} object.
#' @return A data frame of prior parameters for fixed effects.
#' @examples
#' meta <- read.csv(system.file("extdata", "example_metadata.csv.gz", 
#' package = "strainspy"))
#' meta$Case_status <- factor(meta$Case_status)
#' se <- read_sylph(system.file("extdata", "example_sylph_profile.tsv.gz", 
#' package = "strainspy"), meta_data = meta)
#' pri <- suppressWarnings(define_priors(se[1:10, ], 
#' ~ Case_status + Age_at_collection, method = "preset_weak"))
#' head(extract_fixef_priors(pri))
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