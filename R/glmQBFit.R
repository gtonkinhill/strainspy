#' glmQB
#'
#' This function fits a quasibinomial regression model using `stats::glm()`.
#'  
#' It takes a `SummarizedExperiment` object as input, along with a user-defined 
#' formula containing only fixed effect parameters, and fits a quasibinomial 
#' regression model on the assay data.
#'
#' @param se SummarizedExperiment. A `SummarizedExperiment` object containing 
#' the assay data and metadata.
#' @param design Formula. A formula to specify fixed effects, e.g., 
#' ` ~ Age + BMI`.
#' @param nthreads An integer specifying the number of (CPUs or workers) to use. 
#' Defaults to one 1.
#' @param scale_continous Logical. If `TRUE`, all numeric columns in 
#' `colData(se)` are z-score standardized (mean = 0, SD = 1). Default `FALSE`.
#' @param BPPARAM Optional `BiocParallelParam` object. If not provided, the 
#' function will configure an appropriate backend automatically.
#'
#' @return A `strainspy_fit` object with the following components:
#' \item{row_data}{A DFrame with 6 slots with feature details}
#' \item{coefficients}{A DFrame with coefficients for each feature}
#' \item{std_errors}{A DFrame of standard errors for each feature}
#' \item{p_values}{A DFrame of p-values for each feature}
#' \item{residuals}{A DFrame of residual vectors for each feature}
#' \item{convergence}{A named logical vector indicating convergence for each feature}
#' \item{design}{Formula used in the call to `glmQBFit`}
#' \item{call}{Call to `glmQBFit`} 
#'
#' @import SummarizedExperiment
#' 
#'
#' @examples
#' \dontrun{
#' library(strainspy)
#'
#' example_meta_path <- system.file("extdata", "example_metadata.csv.gz", package = "strainspy")
#' example_meta <- readr::read_csv(example_meta_path)
#' example_path <- system.file("extdata", "example_sylph_profile.tsv.gz", package = "strainspy")
#' se <- read_sylph(example_path, example_meta)
#' se <- filter_by_presence(se)
#'
#' design <- as.formula(" ~ Case_status + Age_at_collection")
#'
#' fit <- glmQB(se,  design, nthreads=parallel::detectCores())
#' top_hits(fit, alpha=0.5)
#' plot_manhattan(fit)
#'
#' }
#'
#' @export
glmQBFit <- function(se, design, nthreads=1L, scale_continous=TRUE, BPPARAM=NULL) {
  
  # Validate input
  if (!inherits(se, "SummarizedExperiment")) {
    stop("`se` must be a SummarizedExperiment object.")
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
  
  # check if formula is valid
  nbd = nobars_(design)
  if(is.null(nbd)){
    stop(paste(paste(design, collapse = ''), "--- is not a valid formula."))
  } 
  
  if (!isTRUE(all.equal(nbd, design))) {
    stop(
      paste0(
        "The model formula `", paste(design, collapse = ""), "` contains random effect terms, which are unsupported in `glmQBFit()`. \nConsider using:\n",
        "  - `glmObFit()`\n",
        "  - `glmZiBFit()`\n"
      )
    )
  }
  
  combined_formula <- as.formula(paste(c("Value", as.character(design)),
                                       collapse = " "))
  
  # Ensure that rownames of colData match colnames of the assay
  if (!all(colnames(SummarizedExperiment::assays(se)[[1]]) %in% rownames(col_data))) {
    stop("Column names of assay data do not match row names of colData.")
  }
  
  
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
    function(row_indices) fit_qb_model(SummarizedExperiment::assay(se)[row_indices,],
                                       col_data, combined_formula),
    BPPARAM = BPPARAM
  )
  
  # Flatten the results by removing the first layer of lists
  results <- unname(do.call(c, results))
  
  # Clean up
  BiocParallel::bpstop(BPPARAM)
  
  # sometimes results can be an empty list, remove those dynamically
  rmidx = which(sapply(results, length) == 0)
  
  if(length(rmidx) > 0) {
    seRD = SummarizedExperiment::rowData(se)[-rmidx, , drop = FALSE]
    results[rmidx] <- NULL
  } else {
    seRD = SummarizedExperiment::rowData(se)
  }
  
  # Create the strainspy_fit object
  QBGLM <- new("strainspy_fit",
               row_data = seRD,
               coefficients = DataFrame(purrr::map_dfr(results, ~ .x[[1]][,1])),
               std_errors = DataFrame(purrr::map_dfr(results, ~ .x[[1]][,2])),
               p_values = DataFrame(purrr::map_dfr(results, ~ .x[[1]][,4])),
               zi_coefficients = NULL,
               zi_std_errors = NULL,
               zi_p_values = NULL,
               residuals = DataFrame(purrr::map_dfr(results, ~ .x[[3]])),
               convergence = purrr::map_lgl(results, ~ .x$convergence),
               design = design,
               # assay = assays(se)[[1]],  # Retrieve assay data matrix from SummarizedExperiment
               call = match.call()  # Store the function call for reproducibility
  )
  
  return(QBGLM)
}


#' Fit quasibinomial regression for a Single Feature
#'
#' Fits a quasibinomial regression for a single feature in an assay
#' using `stats::glm()`.
#'
#' @param se_subset A `SummarizedExperiment` object containing the assay data.
#' @param col_data A data frame containing the design matrix and additional covariates.
#' @param combined_formula The formula for the conditional mean model.
#' 
#' #' @return A list with model coefficients, zero-inflation coefficients, residuals,
#'         log-likelihood, and convergence status.
fit_qb_model <- function(se_subset, col_data, combined_formula) {
  
  chunk_results <- lapply(seq_len(nrow(se_subset)), function(row_index){
    # Extract the values for the current feature
    col_data$Value <- base::pmin(as.vector(se_subset[row_index, ]) / 100, 0.99999)
    
    # Run the zero-inflated beta regression
    fit <- tryCatch({
      stats::glm(
        formula = combined_formula,
        data = col_data,
        family = 'quasibinomial'
      )
    }, error = function(e) NULL)
    
    # Handle the case where the model could not be fitted
    if (is.null(fit)) {
      warning("Failed to fit the model for species index: ", row_index)
      return(NULL)
    }
    
    # Extract summary statistics
    smry <- summary(fit)
    
    # Return results as a list
    return(list(
      coefficients = smry$coefficients,
      coefficients_zi = NULL,
      residuals = stats::residuals(fit, type = "response"),
      # log_likelihood = fit$fit$objective,
      convergence = fit$converged
    ))
  })
  
  return(chunk_results)
}


