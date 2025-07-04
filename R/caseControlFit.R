#' caseControlFit
#'
#' This function fits a logistic regression model similar to that published in the original Sylph paper.
#' It takes a `SummarizedExperiment` object as input, along with a user-defined formula.
#'
#' @param se SummarizedExperiment. A `SummarizedExperiment` object containing the assay data and metadata.
#' @param design Formula. A formula to specify the fixed and random effects, e.g., `Outcome ~ Value + Covariate + (1|Sample)`.
#' @param min_identity A numeric value specifying the minimum identity threshold to consider (default=0.98).
#' @param nthreads An integer specifying the number of (CPUs or workers) to use. Defaults
#'        to one 1.
#' @param scale_continous Logical. If `TRUE`, all numeric columns in `colData(se)` are z-score standardized (mean = 0, SD = 1). Defaults to `FALSE`.
#' @param BPPARAM Optional `BiocParallelParam` object. If not provided, the function
#'        will configure an appropriate backend automatically.
#'
#' @return A `strainspy_fit` object with the following components:
#' \item{row_data}{A DFrame with 6 slots with feature details}
#' \item{coefficients}{A DFrame with coefficients for each feature}
#' \item{std_errors}{A DFrame of standard errors for each feature}
#' \item{p_values}{A DFrame of p-values for each feature}
#' \item{residuals}{A DFrame of residual vectors for each feature}
#' \item{convergence}{A named logical vector indicating convergence for each feature}
#' \item{design}{Formula used in the call to `caseControlFit`}
#' \item{call}{Call to `caseControlFit`} 
#'
#' @import SummarizedExperiment
#' @importFrom glmmTMB glmmTMB
#' @importFrom stats terms
#'
#' @examples
#' \dontrun{
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
#' fit <- caseControlFit(se, design, nthreads=parallel::detectCores())
#' min(fit@p_values[,2], na.rm=TRUE)
#'
#' }
#'
#' @export
caseControlFit <- function(se, design, min_identity=0.98, nthreads=1, scale_continous=TRUE, BPPARAM=NULL) {

  # Validate input
  if (!inherits(se, "SummarizedExperiment")) {
    stop("`se` must be a SummarizedExperiment object.")
  }
  
  if(! "Value" %in% attr(stats::terms(design), "term.labels")) {
    stop("`design` must contain the term Value on the RHS, e.g. `Outcome ~ Value + Covariate + (1|Sample)`")
  }

  # colData (sample metadata)
  col_data <- SummarizedExperiment::colData(se)
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
  if (!all(colnames(SummarizedExperiment::assays(se)[[1]]) %in% rownames(col_data))) {
    stop("Column names of assay data do not match row names of colData.")
  }

  # Define priors
  col_data$Value <- 1 # dummy variable for the response
  
  nbd = nobars_(design)
  if(is.null(nbd)){
    stop(paste(paste(design, collapse = ''), "--- is not a valid formula."))
  } 
  nbeta <- ncol(model.matrix(nbd, col_data))
  
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
  results <- unname(do.call(c, results))

  # Clean up
  BiocParallel::bpstop(BPPARAM)

  # Create the strainspy_fit object
  CCGLM <- new("strainspy_fit",
                   row_data = SummarizedExperiment::rowData(se),
                   coefficients = S4Vectors::DataFrame(purrr::map_dfr(results, ~ .x[[1]][,1])),
                   std_errors = S4Vectors::DataFrame(purrr::map_dfr(results, ~ .x[[1]][,2])),
                   p_values = S4Vectors::DataFrame(purrr::map_dfr(results, ~ .x[[1]][,4])),
                   zi_coefficients = NULL,
                   zi_std_errors = NULL,
                   zi_p_values = NULL,
                   residuals = NULL,
                   convergence = purrr::map_lgl(results, ~ .x$convergence),
                   design = design,
                   # assay = assays(se)[[1]],  # Retrieve assay data matrix from SummarizedExperiment
                   call = match.call()  # Store the function call for reproducibility
  )

  return(CCGLM)
}


#' Fit Logistic Regression for a Single Feature
#'
#' Fits a logistic regression for a single feature in an assay
#' using `glmmTMB`.
#'
#' @param se_subset A `SummarizedExperiment` object containing the assay data.
#' @param col_data A data frame containing the design matrix and additional covariates.
#' @param design Formula. A formula to specify the fixed and random effects, e.g., ` ~ Group + (1|Sample)`.
#' @param fixed_priors Optional priors for the model.
#' @param min_identity A numeric value specifying the minimum identity threshold to consider (default=0.98).
#' @return A list with model coefficients, zero-inflation coefficients, residuals,
#'         log-likelihood, and convergence status.
fit_logit_model <- function(se_subset, col_data, design, fixed_priors, min_identity) {

  chunk_results <- lapply(seq_len(nrow(se_subset)), function(row_index){
    # Extract the values for the current feature
    col_data$Value <- base::pmin(as.vector(se_subset[row_index, ]) / 100, 0.99999)
    col_data$Value[col_data$Value < min_identity] <- min_identity
    # Scale to be between 0 and 1. (min-max scaling)
    col_data$Value <- (col_data$Value-min_identity)/(1-min_identity)

    # Run the case-control fit model
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
      convergence = fit$fit$convergence==0
    ))
  })

  return(chunk_results)
}


