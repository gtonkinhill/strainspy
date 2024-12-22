#' glmZiBFit
#'
#' This function fits a zero-inflated beta regression model using the `glmmTMB` package.
#' It takes a `SummarizedExperiment` object as input, along with a user-defined formula,
#' and fits a zero-inflated beta regression model on the assay data.
#'
#' @param se SummarizedExperiment. A `SummarizedExperiment` object containing the assay data and metadata.
#' @param formula Formula. A formula to specify the fixed and random effects, e.g., ` ~ Group + (1|Sample)`.
#' @param assay_name Character. The name of the assay to use from the SummarizedExperiment.
#' Must be one of `"Adjusted_ANI"` or `"Taxonomic_abundance"`.
#' Default is `"Adjusted_ANI"`.
#'
#' @return A list with the following components:
#' \item{coefficients}{A data frame with coefficients, standard errors, z-values, p-values, and FDR for each feature.}
#' \item{fit}{The fitted `glmmTMB` model object.}
#' \item{call}{The original function call.}
#'
#' @import SummarizedExperiment
#' @importFrom glmmTMB glmmTMB
#'
#' @examples
#' \dontrun{
#' library(SummarizedExperiment)
#' library(glmmTMB)
#'
#' example_meta_path <- system.file("extdata", "example_metadata.csv.gz", package = "strainseekr")
#' example_meta <- readr::read_csv(example_meta_path)
#' example_path <- system.file("extdata", "example_sylph_profile.tsv.gz", package = "strainseekr")
#' se <- read_sylph(example_path, example_meta)
#' se <- filter_by_presence(se)
#'
#' design <- as.formula(" ~ Case_status + Age_at_collection")
#' assay_name <- "Adjusted_ANI"
#'
#' results <- glmZiBFit(se,  design)
#' results$coefficients
#' summary(results)
#' }
glmZiBFit <- function(se, design, assay_name = "Adjusted_ANI") {
  # Check if glmmTMB is installed
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    stop("The 'glmmTMB' package is required but is not installed. Please install it with install.packages('glmmTMB').")
  }

  # Validate input
  if (!inherits(se, "SummarizedExperiment")) {
    stop("`se` must be a SummarizedExperiment object.")
  }

  if (!assay_name %in% c("Adjusted_ANI", "Taxonomic_abundance")) {
    stop(paste("The specified assay", assay_name, "is not one of 'Adjusted_ANI' or 'Taxonomic_abundance'."))
  }

  # Extract the assay data (e.g., counts, measurements) and colData (sample metadata)
  assay_data <- assays(se)[[assay_name]]
  col_data <- as.data.frame(colData(se))

  # Ensure the design is a formula
  if (!inherits(design, "formula")) {
    stop("`design` must be a formula (e.g., ~ batch + condition).")
  }
  combined_formula <- as.formula(paste("Value", deparse(design), sep = " "))

  # Ensure that rownames of colData match colnames of the assay
  if (!all(colnames(assay_data) %in% rownames(col_data))) {
    stop("Column names of assay data do not match row names of colData.")
  }

  # Define priors
  nbeta <- ncol(model.matrix(design, col_data))
  fixed_priors <- data.frame(
    prior = rep("normal(0,5)", 2*nbeta),
    class = rep(c("fixef", "fixef_zi"), each=nbeta),
    coef  = rep(as.character(seq(1,nbeta)), 2))

  static_fit <- NULL

  results <- pbapply::pblapply(rownames(assay_data), function(feature) {
    # Extract the values for the current feature (one row of the assay)
    col_data$Value <- base::pmin(as.vector(assay_data[feature, ])/100,
                                 0.99999)

    # Run the zero-inflated beta regression for the current feature
    fit <- tryCatch({
      # Use static_fit for faster convergence, if available
      if (!is.null(static_fit)) {
        update(static_fit, data = col_data)
      } else {
        glmmTMB::glmmTMB(
          formula = combined_formula,
          ziformula = design,
          data = col_data,
          priors = fixed_priors,
          family = glmmTMB::beta_family(link = "logit")
        )
      }
    }, error = function(e) NULL)

    # Retry without static_fit if the first attempt failed
    if (is.null(fit)) {
      print("Refitting without using previous fit!")
      static_fit <<- NULL  # Reset static_fit

      fit <- tryCatch({
        glmmTMB::glmmTMB(
          formula = combined_formula,
          ziformula = design,
          data = col_data,
          priors = fixed_priors,
          family = glmmTMB::beta_family(link = "logit")
        )
      }, error = function(e) NULL)
    }

    # Handle the case where the model could not be fitted at all
    if (is.null(fit)) {
      warning("Failed to fit the model for feature: ", feature)
      return(NULL)
    }

    # Update static_fit with the successful fit for reuse
    if (!is.null(fit)) {
      static_fit <<- fit
    }


    smry <- summary(fit)

    list(
      coefficients <- smry$coefficients$cond,
      coefficients_zi <- smry$coefficients$zi,
      df.residual <- residuals(fit, type = "response"),
      llk <- fit$fit$objective,
      convergence <- fit$fit$convergence
    )
  })

  # Create the betaGLM object
  ZIBetaGLM <- new("betaGLM",
                   coefficients = DataFrame(
                     purrr::map_dfr(results, ~ .x[[1]][,1]) %>%
                       tibble::add_column(.before = 1, strain = rownames(assay_data))
                   ),
                   std_errors = DataFrame(
                     purrr::map_dfr(results, ~ .x[[1]][,2]) %>%
                       tibble::add_column(.before = 1, strain = rownames(assay_data))
                   ),
                   p_values = DataFrame(
                     purrr::map_dfr(results, ~ .x[[1]][,4]) %>%
                       tibble::add_column(.before = 1, strain = rownames(assay_data))
                   ),
                   zi_coefficients = DataFrame(
                     purrr::map_dfr(results, ~ .x[[2]][,1]) %>%
                       tibble::add_column(.before = 1, strain = rownames(assay_data))
                   ),
                   zi_std_errors = DataFrame(
                     purrr::map_dfr(results, ~ .x[[2]][,2]) %>%
                       tibble::add_column(.before = 1, strain = rownames(assay_data))
                   ),
                   zi_p_values = DataFrame(
                     purrr::map_dfr(results, ~ .x[[2]][,4]) %>%
                       tibble::add_column(.before = 1, strain = rownames(assay_data))
                   ),
                   # fdr = NULL,  # Placeholder for FDR values, update as needed
                   # zi_fdr = NULL,  # Placeholder for Zero-Inflated FDR values, update as needed
                   residuals = DataFrame(
                     purrr::map_dfr(results, ~ .x[[3]]) %>%
                       tibble::add_column(.before = 1, strain = rownames(assay_data))
                   ),
                   design = model.matrix(design, data = as.data.frame(colData(se))),
                   assay = assays(se)[[assay_name]],  # Retrieve assay data matrix from SummarizedExperiment
                   call = match.call()  # Store the function call for reproducibility
  )

  return(ZIBetaGLM)
}
