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
#' @param scale_continous Logical. If `TRUE`, all numeric columns in `colData(se)` are z-score standardized (mean = 0, SD = 1). Defaults to `FALSE`.
#' @param BPPARAM Optional `BiocParallelParam` object. If not provided, the function
#'        will configure an appropriate backend automatically.
#' @param method Character. The method to use for fitting the model. Either 'glmmTMB' (default) or 'gamlss'.
#'
#' @return A list with the following components:
#' \item{coefficients}{A data frame with coefficients, standard errors, z-values, p-values, and FDR for each feature.}
#' \item{fit}{The fitted `glmmTMB` model object.}
#' \item{call}{The original function call.}
#'
#' @import SummarizedExperiment
#' @import gamlss
#' @importFrom glmmTMB glmmTMB
#' @importFrom stats model.matrix
#' @importFrom methods new
#' @importFrom stats as.formula
#' 
#'
#' @examples
#' \dontrun{
#' library(SummarizedExperiment)
#' library(strainspy)
#' library(glmmTMB)
#' library(gamlss)
#'
#' example_meta_path <- system.file("extdata", "example_metadata.csv.gz", package = "strainspy")
#' example_meta <- readr::read_csv(example_meta_path)
#' example_path <- system.file("extdata", "example_sylph_profile.tsv.gz", package = "strainspy")
#' se <- read_sylph(example_path, example_meta)
#' se <- filter_by_presence(se)
#'
#' design <- as.formula(" ~ Case_status + Age_at_collection + (1|Sex)")
#' design <- as.formula(" ~ Case_status + Age_at_collection + gamlss::random(Sex)")
#' design <- as.formula(" ~ Case_status")
#'
#' fit <- glmZiBFit(se[1:5,], design, nthreads=1, method='gamlss')
#' top_hits(fit, alpha=1)
#'
#' }
#'
#' @export
glmZiBFit <- function(se, design, nthreads=1, scale_continous=TRUE, BPPARAM=NULL,
                      method='glmmTMB') {
  # Check if glmmTMB is installed
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    stop("The 'glmmTMB' package is required but is not installed. Please install it with install.packages('glmmTMB').")
  }

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

  # Define priors
  nbeta <- ncol(model.matrix(nbd, col_data))
  fixed_priors <- data.frame(
    prior = rep("normal(0,5)", 2*nbeta),
    class = rep(c("fixef", "fixef_zi"), each=nbeta),
    coef  = rep(as.character(seq(1,nbeta)), 2))

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

  if (method=='glmmTMB'){
    results <- BiocParallel::bplapply(
      row_chunks,
      function(row_indices) fit_zero_inflated_beta(SummarizedExperiment::assay(se)[row_indices, , drop=FALSE],
                                                   col_data, combined_formula, design, fixed_priors),
      BPPARAM = BPPARAM
    )
  } else{
    results <- BiocParallel::bplapply(
      row_chunks,
      function(row_indices) fit_zero_inflated_beta_gamlss(SummarizedExperiment::assay(se)[row_indices, , drop=FALSE],
                                                   col_data, combined_formula, design, fixed_priors),
      BPPARAM = BPPARAM
    )
  }


  # Flatten the results by removing the first layer of lists
  results <- unname(do.call(c, results))

  # Clean up
  BiocParallel::bpstop(BPPARAM)
  
  # sometimes results can be an empty list, remove those dynamically
  rmidx = which(sapply(results, length) == 0)
  # sometimes the loglikelihood is NA when the result failed
  rmidx = c(rmidx, which(is.na(unlist(lapply(results, function(x) x$log_likelihood)))))
  rmidx = unique(rmidx)
  if(length(rmidx) > 0) {
    seRD = SummarizedExperiment::rowData(se)[-rmidx, , drop = FALSE]
    results[rmidx] <- NULL
  } else {
    seRD = SummarizedExperiment::rowData(se)
  }

  # Create the betaGLM object
  ZIBetaGLM <- methods::new("betaGLM",
                            row_data = seRD,
                            coefficients = DataFrame(purrr::map_dfr(results, ~ .x[[1]][,1])),
                            std_errors = DataFrame(purrr::map_dfr(results, ~ .x[[1]][,2])),
                            p_values = DataFrame(purrr::map_dfr(results, ~ .x[[1]][,4])),
                            zi_coefficients = DataFrame(purrr::map_dfr(results, ~ .x[[2]][,1])),
                            zi_std_errors = DataFrame(purrr::map_dfr(results, ~ .x[[2]][,2])),
                            zi_p_values = DataFrame(purrr::map_dfr(results, ~ .x[[2]][,4])),
                            residuals = DataFrame(purrr::map_dfr(results, ~ .x[[3]])),
                            convergence = purrr::map_lgl(results, ~ .x$convergence),
                            design = design,
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
    col_data$Value <- base::pmin(as.vector(se_subset[row_index, ]) / 100, 0.99999)

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
      convergence = fit$fit$convergence==0
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
#' @return A list with model coefficients, zero-inflation coefficients, residuals,
#'         log-likelihood, and convergence status.
#'
#' @importFrom stats residuals
fit_zero_inflated_beta_gamlss <- function(se_subset, col_data, combined_formula, design, fixed_priors) {

  chunk_results <- lapply(seq_len(nrow(se_subset)), function(row_index){
    # Extract the values for the current feature
    col_data$Value <- base::pmin(as.vector(se_subset[row_index, ]) / 100, 0.99999)
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
      rep(.x, sum(!grepl('random', names(fit[[paste0(.x, '.coefficients')]]))))
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
      residuals = NULL,
      log_likelihood = NA,
      convergence = fit$converged
    ))
  })

  return(chunk_results)
}


