#' abundanceFit
#'
#' This function fits a linear regression model to proportional abundance data using either `stats::lm()` or `lmerTest::lmer()` for models with random effects. 
#' It takes a `SummarizedExperiment` object as input, along with a user-defined formula,
#' and fits a linear model on the assay data. It transforms the abundance data using: `log10(proportional_abundance + 1e-5)` prior to fitting.
#'
#' @param se SummarizedExperiment. A `SummarizedExperiment` object containing the assay data and metadata.
#' @param design Formula. A formula to specify the fixed and random effects, e.g., ` ~ Group + (1|Sample)`.
#' @param nthreads An integer specifying the number of (CPUs or workers) to use. Defaults
#'        to one 1.
#' @param scale_continous Logical. If `TRUE`, all numeric columns in `colData(se)` are z-score standardized (mean = 0, SD = 1). Defaults to `FALSE`.
#' @param transform If data is already transformed, set to `NULL` (default). Supported options: arcsin transform `arcsin` or centered log ratio `CLR`. `CLR` requires `compositions` package.
#' @param BPPARAM Optional `BiocParallelParam` object. If not provided, the function
#'        will configure an appropriate backend automatically.
#'
#' @return A list with the following components:
#' \item{coefficients}{A data frame with coefficients, standard errors, z-values, p-values, and FDR for each feature.}
#' \item{fit}{The fitted `glmmTMB` model object.}
#' \item{call}{The original function call.}
#' 
#' @importFrom lmerTest lmer
#' @import SummarizedExperiment
#' @importFrom glmmTMB glmmTMB
#' @importFrom stats residuals logLik
#'
#' @examples
#' \dontrun{
#' library(SummarizedExperiment)
#' library(strainspy)
#' library(glmmTMB)
#'
#' example_meta_path <- system.file("extdata", "example_metadata.csv.gz", package = "strainspy")
#' example_meta <- readr::read_csv(example_meta_path)
#' example_path <- system.file("extdata", "example_sylph_profile.tsv.gz", package = "strainspy")
#' se <- read_sylph(example_path, example_meta)
#' se <- filter_by_presence(se)
#'
#' design <- as.formula(" ~ Case_status + Age_at_collection")
#'
#' fit <- glmFit(se,  design, nthreads=4, family=glmmTMB::ordbeta())
#' summary(fit)
#'
#' }
#'
#' @export
abundanceFit <- function(se, design, nthreads=1, scale_continous=TRUE, transform=NULL, BPPARAM=NULL) {
  # Check if glmmTMB is installed
  if (!requireNamespace("lmerTest", quietly = TRUE)) {
    stop("The 'lmerTest' package is required but is not installed. Please install it with install.packages('lmerTest').")
  }
  
  # Validate input
  if (!inherits(se, "SummarizedExperiment")) {
    stop("`se` must be a SummarizedExperiment object.")
  }
  
  # colData (sample metadata)
  col_data <- as.data.frame(SummarizedExperiment::colData(se))
  if (scale_continous==TRUE){
    for (col in names(col_data)) {
      if (is.numeric(col_data[[col]])) {
        col_data[[col]] <- as.numeric(scale(col_data[[col]]))  # Scale numeric columnsa + convert from matrix to vector (we are guaranteed to go col by col)
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
  
  # transform - this is probably not the optimal way to do this!
  if(!is.null(transform)){
    if(transform == "arcsin") {
      
      SummarizedExperiment::assay(se) <- as(apply(SummarizedExperiment::assay(se) , 2, function(x) asin(sqrt(x/100))), 'dgRMatrix')
      
    } else if(transform == "CLR") {
      if (!requireNamespace("compositions", quietly = TRUE)) {
        stop("The 'compositions' package is required for CLR transformation, but is not installed. Please install it with install.packages('compositions').")
      }
      
      SummarizedExperiment::assay(se) <- as(apply(SummarizedExperiment::assay(se) , 2, function(x) compositions::clr(x)), 'dgRMatrix')
      
    } else {
      stop("Unknown transform = ", transform, ' --- see ?strainspy::abundanceFit')
    }
  }
  
  # simple model
  if( length(grep("\\|", deparse(design))) == 0 ) {
    results <- BiocParallel::bplapply(
      row_chunks,
      function(row_indices) fit_lm(se_subset = SummarizedExperiment::assay(se)[row_indices,],
                                   col_data = col_data, 
                                   combined_formula =  combined_formula),
      BPPARAM = BPPARAM
    )
  } else {
    results <- BiocParallel::bplapply(
      row_chunks,
      function(row_indices) fit_lmer(se_subset = SummarizedExperiment::assay(se)[row_indices,],
                                     col_data = col_data, 
                                     combined_formula =  combined_formula),
      BPPARAM = BPPARAM
    )
  }
  
  
  
  # Flatten the results by removing the first layer of lists
  results <- do.call(c, results)
  
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
  
  # Create the betaGLM object
  p_val_idx = ncol(results[[1]]$coefficients)
  
  GLM <- new("betaGLM",
             row_data = seRD,
             coefficients = DataFrame(purrr::map_dfr(results, ~ .x[[1]][,1])),
             std_errors = DataFrame(purrr::map_dfr(results, ~ .x[[1]][,2])),
             p_values = DataFrame(purrr::map_dfr(results, ~ .x[[1]][,p_val_idx])),
             zi_coefficients = NULL,
             zi_std_errors = NULL,
             zi_p_values = NULL,
             residuals = DataFrame(purrr::map_dfr(results, ~ .x[[3]])),
             convergence = purrr::map_lgl(results, ~ .x$convergence),
             design = design,
             # assay = assays(se)[[1]],  # Retrieve assay data matrix from SummarizedExperiment
             call = match.call()  # Store the function call for reproducibility
  )
  
  return(GLM)
}



#' Fit a linear model for a Single Feature
#'
#' Fits a simple linear model for functions without random effects using stats::lm()
#'
#' @param se A `SummarizedExperiment` object containing the assay data.
#' @param col_data A data frame containing the design matrix and additional covariates.
#' @param combined_formula The formula for the conditional mean model.
#' @return A list with model coefficients, zero-inflation coefficients, residuals,
#'         log-likelihood, and convergence status.
fit_lm <- function(se_subset, col_data, combined_formula) {
  
  chunk_results <- lapply(seq_len(nrow(se_subset)), function(row_index){
    # Extract the values for the current feature
    col_data$Value <- se_subset[row_index, ] 
    # for ordBeta, values can be in the closed interval [0,1], but for log linked functions, best to avoid 0
    
    # Run the zero-inflated beta regression
    fit <- tryCatch({
      stats::lm(
        formula = combined_formula,
        data = col_data
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
      log_likelihood = as.numeric(stats::logLik(fit)),
      # OLS usually converged when solution exists.
      convergence = 0
    ))
  })
  
  return(chunk_results)
}

#' Fit a linear model with random effects for a single feature
#'
#' Fits a linear model for a single feature in an assay using `lmerTest::lmer()`
#'
#' @param se A `SummarizedExperiment` object containing the assay data.
#' @param col_data A data frame containing the design matrix and additional covariates.
#' @param combined_formula The formula for the conditional mean model.
#' @return A list with model coefficients, zero-inflation coefficients, residuals,
#'         log-likelihood, and convergence status.
fit_lmer <- function(se_subset, col_data, combined_formula) {
  
  chunk_results <- lapply(seq_len(nrow(se_subset)), function(row_index){
    # Extract the values for the current feature
    col_data$Value <- se_subset[row_index, ]
    # for ordBeta, values can be in the closed interval [0,1], but for log linked functions, best to avoid 0
    
    # Run the zero-inflated beta regression
    fit <- tryCatch({
      suppressMessages(
        suppressWarnings(
          lmerTest::lmer(
            formula = combined_formula,
            data = col_data
          )
        )
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
      log_likelihood = as.numeric(stats::logLik(fit)),
      # OLS usually converged when solution exists.
      convergence = fit@optinfo$conv$opt
    ))
  })
  
  return(chunk_results)
}


