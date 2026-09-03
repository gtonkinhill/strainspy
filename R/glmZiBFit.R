#' glmZiBFit
#'
#' This function fits a zero-inflated beta regression model using either `glmmTMB` or `gamlss` packages.
#' It takes a `SummarizedExperiment` object as input, along with a user-defined formula,
#' and fits a zero-inflated beta regression model on the assay data.
#'
#' @param se SummarizedExperiment. A `SummarizedExperiment` object containing the assay data and metadata.
#' @param design Formula. A formula to specify the fixed and random effects, e.g., ` ~ Group + (1|Sample)`.
#' @param nthreads An integer specifying the number of (CPUs or workers) to use. Defaults
#'        to one 1.
#' @param scale_continuous Logical. If `TRUE`, all numeric columns in `colData(se)` are z-score standardized (mean = 0, SD = 1). Defaults to `TRUE`.
#' @param MAP_prior One of `strainspy_priors` object, character string (`"preset_weak", "preset_strong"`),
#' a `data.frame` of priors, or `NULL`. Default `"preset_weak"`. See Details. 
#' @param BPPARAM Optional `BiocParallelParam` object. If not provided, the function
#'        will configure an appropriate backend automatically.
#' @param method Character. The method to use for fitting the model. Either 'glmmTMB' (default) or 'gamlss'.
#' @param return_vcov Logical. Return variance-covariance matrices for each model, required for `estimate_effect_sizes()`. 
#' Disabling will result in a much smaller `strainspy_fit` object that is useful for most analyses. Defaults to `TRUE`.
#' @param progress Logical. If `TRUE` (default), progress bars and progress messages are printed. Set to `FALSE` to silence them.
#'
#' @return A `strainspy_fit` object with the following components:
#' \item{row_data}{A DFrame with 6 slots with feature details}
#' \item{coefficients}{A DFrame with coefficients for each feature}
#' \item{std_errors}{A DFrame of standard errors for each feature}
#' \item{p_values}{A DFrame of p-values for each feature}
#' \item{zi_coefficients}{A DFrame of zero inflated coefficients for each feature}
#' \item{zi_std_errors}{A DFrame of zero inflated standard errors for each feature}
#' \item{zi_p_values}{A DFrame of zero inflated p-values for each feature}
#' \item{residuals}{A DFrame of residual vectors for each feature}
#' \item{convergence}{A named logical vector indicating convergence for each feature}
#' \item{design}{Formula used in the call to `glmZiBFit`}
#' \item{call}{Call to `glmZiBFit`} 
#'
#' @import SummarizedExperiment
#' @importFrom glmmTMB glmmTMB
#' @importFrom stats model.matrix
#' @importFrom methods new
#' @importFrom stats as.formula
#' 
#'
#' @examples
#' library(strainspy)
#'
#' example_meta_path <- system.file("extdata", "example_metadata.csv.gz", 
#' package = "strainspy")
#' example_meta <- read.csv(example_meta_path)
#' example_path <- system.file("extdata", "example_sylph_profile.tsv.gz", 
#' package = "strainspy")
#' se <- read_sylph(example_path, example_meta)
#' se <- filter_by_presence(se)
#'
#' design <- as.formula(" ~ Case_status")
#'
#' fit <- glmZiBFit(se[1:10,], design, nthreads=2)
#'
#' @export
glmZiBFit <- function(se, design, nthreads=1, scale_continuous=TRUE, 
                      MAP_prior = 'preset_weak', BPPARAM=NULL, method='glmmTMB',
                      return_vcov = TRUE, progress=TRUE) {
  check_progress(progress)
  method <- match.arg(method, choices = c("glmmTMB", "gamlss"))

  # Check if glmmTMB is installed when selected as backend
  if (method == "glmmTMB" && !requireNamespace("glmmTMB", quietly = TRUE)) {
    stop("The 'glmmTMB' package is required for method = 'glmmTMB' but is not installed. Please install it with install.packages('glmmTMB').")
  }

  # Check if gamlss is installed when selected as backend
  if (method == "gamlss") {
    missing_pkgs <- c("gamlss", "gamlss.dist")[
      !vapply(c("gamlss", "gamlss.dist"), requireNamespace, logical(1), quietly = TRUE)
    ]
    if (length(missing_pkgs) > 0) {
      stop("The following package(s) are required for method = 'gamlss' but are not installed: ",
           paste(missing_pkgs, collapse = ", "),
           ". Please install them with install.packages(c(",
           paste0('"', missing_pkgs, '"', collapse = ", "), ")).")
    }
  }
  
  # Validate input
  if (!inherits(se, "SummarizedExperiment")) {
    stop("`se` must be a SummarizedExperiment object.")
  }
  
  # colData (sample metadata)
  col_data <- colData(se)
  if (scale_continuous==TRUE){
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
  
  combined_formula <- as.formula(paste(c("Value", as.character(design)),
                                       collapse = " "))
  
  # Ensure that rownames of colData match colnames of the assay
  if (!all(colnames(assays(se)[[1]]) %in% rownames(col_data))) {
    stop("Column names of assay data do not match row names of colData.")
  }
  
  prior_obj = resolve_priors(MAP_prior, se, nbd)
  # prior_obj is NULL when MAP_prior = NULL (unpenalised / MLE fit)
  fixed_priors = if (is.null(prior_obj)) NULL else prior_obj@priors_df
  
  # Set up parallel infrastructure
  if ((nthreads > 1) & (.Platform$OS.type != "windows")) {
    # Check the operating system and set the backend accordingly
    if (is.null(BPPARAM)) {
      # Use MulticoreParam for Unix-like systems
      # BPPARAM <- BiocParallel::MulticoreParam(
      #   workers = nthreads
      # )
      BPPARAM <- BiocParallel::SnowParam(workers = nthreads, progressbar = progress, tasks=100)
    } else if (!progress) {
      BiocParallel::bpprogressbar(BPPARAM) <- FALSE
    }
  } else {
    BPPARAM <- BiocParallel::SerialParam(progressbar = progress)
  }
  
  # Split rows into 50-row chunks
  row_chunks <- split(
    seq_len(nrow(se)),
    ceiling(seq_len(nrow(se)) / 100)  # 50 rows per chunk
  )
  
  if (progress) message("Fitting model...")
  
  if (method == 'glmmTMB') {
    results <- BiocParallel::bplapply(
      row_chunks,
      function(row_indices) fit_zero_inflated_beta(SummarizedExperiment::assay(se)[row_indices, , drop=FALSE],
                                                   col_data, combined_formula, design, fixed_priors),
      BPPARAM = BPPARAM
    )
  } else {
    results <- BiocParallel::bplapply(
      row_chunks,
      function(row_indices) fit_zero_inflated_beta_gamlss(SummarizedExperiment::assay(se)[row_indices, , drop=FALSE],
                                                          col_data, combined_formula, design, fixed_priors,
                                                          progress = progress),
      BPPARAM = BPPARAM
    )
  }
  
  
  # Flatten the results by removing the first layer of lists
  results <- unname(do.call(c, results))
  
  # Clean up
  BiocParallel::bpstop(BPPARAM)
  
  # Drop features whose fit failed. See `failed_fit_index()` for why both the
  # zero-length and the NA-log_likelihood protocols have to be handled.
  rmidx = failed_fit_index(results)
  if(length(rmidx) > 0) {
    seRD = SummarizedExperiment::rowData(se)[-rmidx, , drop = FALSE]
    results[rmidx] <- NULL
  } else {
    seRD = SummarizedExperiment::rowData(se)
  }
  
  # Create the strainspy_fit object
  # This needs to be triple checked for compatibility with other models and families
  if(return_vcov){
    ZIBetaGLM <- methods::new("strainspy_fit",
                              row_data = seRD,
                              priors = prior_obj,
                              scale_continuous = scale_continuous,
                              coefficients = DataFrame(purrr::map_dfr(results, ~ .x$coefficients[,1])),
                              std_errors = DataFrame(purrr::map_dfr(results, ~ .x$coefficients[,2])),
                              p_values = DataFrame(purrr::map_dfr(results, ~ .x$coefficients[,4])),
                              vcov = List(purrr::map(results, ~ as.data.frame(.x$vcov))),
                              zi_coefficients = DataFrame(purrr::map_dfr(results, ~ .x$coefficients_zi[,1])),
                              zi_std_errors = DataFrame(purrr::map_dfr(results, ~ .x$coefficients_zi[,2])),
                              zi_p_values = DataFrame(purrr::map_dfr(results, ~ .x$coefficients_zi[,4])),
                              zi_vcov = List(purrr::map(results, ~ as.data.frame(.x$zi_vcov))),
                              residuals = DataFrame(purrr::map_dfr(results, ~ .x$residuals)),
                              convergence = purrr::map_lgl(results, ~ .x$convergence),
                              design = design,
                              call = match.call(),  # Store the function call for reproducibility
                              family = glmmTMB::beta_family(link = "logit")
    )
  } else {
    ZIBetaGLM <- methods::new("strainspy_fit",
                              row_data = seRD,
                              priors = prior_obj,
                              scale_continuous = scale_continuous,
                              coefficients = DataFrame(purrr::map_dfr(results, ~ .x$coefficients[,1])),
                              std_errors = DataFrame(purrr::map_dfr(results, ~ .x$coefficients[,2])),
                              p_values = DataFrame(purrr::map_dfr(results, ~ .x$coefficients[,4])),
                              vcov = NULL,
                              zi_coefficients = DataFrame(purrr::map_dfr(results, ~ .x$coefficients_zi[,1])),
                              zi_std_errors = DataFrame(purrr::map_dfr(results, ~ .x$coefficients_zi[,2])),
                              zi_p_values = DataFrame(purrr::map_dfr(results, ~ .x$coefficients_zi[,4])),
                              zi_vcov = NULL,
                              residuals = DataFrame(purrr::map_dfr(results, ~ .x$residuals)),
                              convergence = purrr::map_lgl(results, ~ .x$convergence),
                              design = design,
                              call = match.call(),  # Store the function call for reproducibility
                              family = glmmTMB::beta_family(link = "logit")
    )
  }
  
  
  
  return(ZIBetaGLM)
}


#' Fit Zero-Inflated Beta Regression for a Single Feature
#'
#' Fits a zero-inflated beta regression for a single feature in an assay
#' using `glmmTMB`.
#'
#' @param se_subset A `SummarizedExperiment` object containing the assay data.
#' @param col_data A data frame containing the design matrix and additional covariates.
#' @param combined_formula The formula for the conditional mean model.
#' @param design The formula for the zero-inflation model.
#' @param fixed_priors Optional priors for the model.
#' @return A list with model coefficients, zero-inflation coefficients, residuals,
#'         log-likelihood, and convergence status.
#'
#' @importFrom stats residuals
fit_zero_inflated_beta <- function(se_subset, col_data, combined_formula, design, fixed_priors) {
  chunk_results <- lapply(seq_len(nrow(se_subset)), function(row_index){
    # Extract the values for the current feature
    col_data$Value <- offset_ANI(as.vector(se_subset[row_index, ])/100)
    
    # Run the zero-inflated beta regression
    fit <- tryCatch({
      glmmTMB::glmmTMB(
        formula = combined_formula,
        ziformula = design,
        data = col_data,
        priors = fixed_priors,
        family = glmmTMB::beta_family(link = "logit")
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
        coefficients_zi = na_matrix,
        residuals = rep(NA, nrow(col_data)),
        log_likelihood = NA,
        vcov = NA,
        zi_vcov = NA,
        convergence = FALSE
      ))
      
    }
    
    # Extract summary statistics
    smry <- summary(fit)
    
    # Return results as a list
    return(list(
      coefficients = smry$coefficients$cond,
      coefficients_zi = smry$coefficients$zi,
      residuals = stats::residuals(fit, type = "response"),
      log_likelihood = fit$fit$objective,
      vcov = smry$vcov$cond,
      zi_vcov = smry$vcov$zi,
      convergence = fit$fit$convergence == 0 & isTRUE(fit$sdr$pdHess)
    ))
  })
  
  return(chunk_results)
}



#' Fit Zero-Inflated Beta Regression for a Single Feature using the gamlss package
#'
#' Fits a zero-inflated beta regression for a single feature in an assay
#' using `gamlss`.
#'
#' @param se_subset A `SummarizedExperiment` object containing the assay data.
#' @param col_data A data frame containing the design matrix and additional covariates.
#' @param combined_formula The formula for the conditional mean model.
#' @param design The formula for the zero-inflation model.
#' @param fixed_priors Optional priors for the model.
#' @param progress Logical. If `TRUE` (default), a progress message is printed.
#' @return A list with model coefficients, zero-inflation coefficients, residuals,
#'         log-likelihood, and convergence status.
#'
#' @importFrom stats residuals
fit_zero_inflated_beta_gamlss <- function(se_subset, col_data, combined_formula, design, fixed_priors,
                                          progress = TRUE) {
  if (progress) message("Fitting model...")
  chunk_results <- lapply(seq_len(nrow(se_subset)), function(row_index){
    # Extract the values for the current feature
    col_data$Value <- offset_ANI(as.vector(se_subset[row_index, ])/100)
    
    temp_dat <- col_data[, all.vars(combined_formula)]
    
    # convert characters to factors
    temp_dat[] <- lapply(temp_dat, function(x) {
      if (is.character(x)) factor(x) else x
    })
    
    # Run the zero-inflated beta regression
    fit <- tryCatch({
      gamlss::gamlss(
        formula = combined_formula,
        sigma.fo = ~1,
        nu.fo = design,
        data = temp_dat,
        family = gamlss.dist::BEZI,
        control = gamlss::gamlss.control(trace=FALSE),
        # mu.start = 0.5,
        # nu.start = 0.3,
        # sigma.start = 1,
        # tau.start = 1
      )
    }, error = function(e) {
      return(NULL)
    })
    
    
    # Handle the case where the model could not be fitted
    m <- model.matrix(design, temp_dat)
    na_matrix <- matrix(NA,
                        nrow = ncol(m),
                        ncol = 4,
                        dimnames = list(colnames(m),
                                        c("Estimate", "Std. Error", "z value", "Pr(>|z|)")))
    
    if (is.null(fit)) {
      warning("Failed to fit the model for species index: ", row_index)
      
      return(list(
        coefficients = na_matrix,
        coefficients_zi = na_matrix,
        residuals = rep(NA, nrow(col_data)),
        log_likelihood = NA,
        convergence = FALSE
      ))
    }
    
    # Extract summary statistics
    stdout <- capture.output(smry <- summary(fit, robust=TRUE))
    
    index <- unlist(purrr::map(fit$parameters, ~{
      rep(.x, sum(!grepl('random', names(fit[[paste0(.x, '.coefficients')]]) )))
    }))
    
    if (length(index) != nrow(smry)) {
      warning("Failed to fit the model for species index: ", row_index)
      
      return(list(
        coefficients = na_matrix,
        coefficients_zi = na_matrix,
        residuals = rep(NA, nrow(col_data)),
        log_likelihood = NA,
        convergence = FALSE
      ))
    }
    
    smry <- purrr::map(fit$parameters, ~{
      smry[.x==index,,drop=FALSE]
    })
    
    names(smry) <- fit$parameters
    
    # Return results as a list
    return(list(
      coefficients = smry$mu,
      coefficients_zi = smry$nu,
      residuals = fit$residuals,
      log_likelihood = as.numeric(logLik(fit)),
      convergence = fit$converged
    ))
  })
  
  return(chunk_results)
}