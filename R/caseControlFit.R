#' caseControlFit
#'
#' This function fits a logistic regression model similar to that published in the original Sylph paper.
#' It takes a `SummarizedExperiment` object as input, along with a user-defined formula.
#'
#' @param se SummarizedExperiment. A `SummarizedExperiment` object containing the assay data and metadata.
#' @param formula Formula. A formula to specify the fixed and random effects, e.g., ` ~ Group + (1|Sample)`.
#' @param min_identity A numeric value specifying the minimum identity threshold to consider (default=0.98).
#' @param nthreads An integer specifying the number of (CPUs or workers) to use. Defaults
#'        to one 1.
#' @param family A `glmmTMB` family object. Defaults to `glmmTMB::ordbeta()`.
#' @param BPPARAM Optional `BiocParallelParam` object. If not provided, the function
#'        will configure an appropriate backend automatically.
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
#' library(strainspy)
#'
#' example_meta_path <- system.file("extdata", "example_metadata.csv.gz", package = "strainspy")
#' example_meta <- readr::read_csv(example_meta_path)
#' example_meta$Case_status <- factor(example_meta$Case_status)
#' example_path <- system.file("extdata", "example_sylph_profile.tsv.gz", package = "strainspy")
#' se <- read_sylph(example_path, example_meta)
#' se <- filter_by_presence(se)
#'
#' design <- as.formula("Case_status ~ Value + Age_at_collection")
#'
#' fit <- caseControlFit(se,  design, nthreads=4)
#' summary(fit)
#'
#' }
#'
#' @export
caseControlFit <- function(se, design, min_identity=0.98, nthreads=1, scale_continous=TRUE, BPPARAM=NULL) {

  # Validate input
  if (!inherits(se, "SummarizedExperiment")) {
    stop("`se` must be a SummarizedExperiment object.")
  }

  # colData (sample metadata)
  col_data <- colData(se)
  if (scale_continous==TRUE){
    for (col in names(col_data)) {
      if (is.numeric(col_data[[col]])) {
        col_data[[col]] <- scale(col_data[[col]])  # Scale numeric columns
      }
    }
  }

  # Ensure the design is a formula
  if (!inherits(design, "formula")) {
    stop("`design` must be a formula (e.g., ~ batch + condition).")
  }

  # Ensure that rownames of colData match colnames of the assay
  if (!all(colnames(assays(se)[[1]]) %in% rownames(col_data))) {
    stop("Column names of assay data do not match row names of colData.")
  }

  # Define priors
  col_data$Value <- 1 # dummy variable for the response
  nbeta <- ncol(model.matrix(design, col_data))
  fixed_priors <- data.frame(
    prior = rep("normal(0,5)", nbeta),
    class = rep("fixef", each=nbeta),
    coef  = as.character(seq(1,nbeta)))

  # Set up parallel infrastructure
  if ((nthreads > 1) & (.Platform$OS.type != "windows")) {
    # Check the operating system and set the backend accordingly
    if (is.null(BPPARAM)) {
      # Use MulticoreParam for Unix-like systems
      # BPPARAM <- BiocParallel::MulticoreParam(
      #   workers = nthreads
      # )
      BPPARAM <- BiocParallel::SnowParam(workers = nthreads, progressbar = TRUE, tasks=100)
    }
  } else {
    BPPARAM <- BiocParallel::SerialParam(progressbar = TRUE)
  }

  # Split rows into 50-row chunks
  row_chunks <- split(
    seq_len(nrow(se)),
    ceiling(seq_len(nrow(se)) / 100)  # 50 rows per chunk
  )

  results <- BiocParallel::bplapply(
    row_chunks,
    function(row_indices) fit_logit_model(SummarizedExperiment::assay(se)[row_indices,],
                                    col_data, design, fixed_priors, min_identity),
    BPPARAM = BPPARAM
  )

  # Flatten the results by removing the first layer of lists
  results <- do.call(c, results)

  # Clean up
  BiocParallel::bpstop(BPPARAM)

  # Create the betaGLM object
  ZIBetaGLM <- new("betaGLM",
                   row_data = rowData(se),
                   coefficients = DataFrame(purrr::map_dfr(results, ~ .x[[1]][,1])),
                   std_errors = DataFrame(purrr::map_dfr(results, ~ .x[[1]][,2])),
                   p_values = DataFrame(purrr::map_dfr(results, ~ .x[[1]][,4])),
                   zi_coefficients = NULL,
                   zi_std_errors = NULL,
                   zi_p_values = NULL,
                   residuals = NULL, # DataFrame(purrr::map_dfr(results, ~ .x[[3]])),
                   # design = model.matrix(design, data = as.data.frame(colData(se))),
                   # assay = assays(se)[[1]],  # Retrieve assay data matrix from SummarizedExperiment
                   call = match.call()  # Store the function call for reproducibility
  )

  return(ZIBetaGLM)
}


#' Fit Zero-Inflated Beta Regression for a Single Feature
#'
#' Fits a zero-inflated beta regression for a single feature in an assay
#' using `glmmTMB`.
#'
#' @param se A `SummarizedExperiment` object containing the assay data.
#' @param row_index The index of the feature (row) to be processed.
#' @param col_data A data frame containing the design matrix and additional covariates.
#' @param combined_formula The formula for the conditional mean model.
#' @param fixed_priors Optional priors for the model.
#' @param nt Number of threads for parallel computation in model fitting.
#' @param feature Optional name of the feature for debugging or error messages.
#' @return A list with model coefficients, zero-inflation coefficients, residuals,
#'         log-likelihood, and convergence status.
fit_logit_model <- function(se_subset, col_data, design, fixed_priors, min_identity) {

  chunk_results <- lapply(seq_len(nrow(se_subset)), function(row_index){
    # Extract the values for the current feature
    col_data$Value <- base::pmin(as.vector(se_subset[row_index, ]) / 100, 0.99999)
    col_data$Value[col_data$Value < min_identity] <- min_identity

    # Run the zero-inflated beta regression
    fit <- tryCatch({
      glmmTMB::glmmTMB(
        formula = design,
        data = col_data,
        priors = fixed_priors,
        family = binomial(link = "logit")
      )
    }, error = function(e) NULL)

    # Handle the case where the model could not be fitted
    if (is.null(fit)) {
      warning("Failed to fit the model for species index: ", row_index)

      m <- model.matrix(design, col_data)
      na_matrix <- matrix(NA,
                          nrow = ncol(m),
                          ncol = 4,
                          dimnames = list(colnames(m),
                                          c("Estimate", "Std. Error", "z value", "Pr(>|z|)")))

      return(list(
        coefficients = na_matrix,
        coefficients_zi = NULL,
        residuals = rep(NA, nrow(col_data)),
        log_likelihood = NA,
        convergence = FALSE
      ))
    }

    # Extract summary statistics
    smry <- summary(fit)

    # Return results as a list
    return(list(
      coefficients = smry$coefficients$cond,
      coefficients_zi = NULL,
      residuals = residuals(fit, type = "response"),
      log_likelihood = fit$fit$objective,
      convergence = fit$fit$convergence
    ))
  })

  return(chunk_results)
}


