#' Constructor for strainspy_priors using empirical Bayes estimation
#'
#' Fits fast binomial and beta models to estimate MAP prior SDs for use in 
#' zero-inflated and ordinal beta models fitted using `GLMMTMB`. 
#' 
#' @importFrom stats quantile
#' 
#' @param se A \code{SummarizedExperiment} object with count matrix and colData.
#' @param design A fixed-effects-only formula (e.g., ~ group + age + sex).
#' @param nthreads An integer specifying the number of (CPUs or workers) to use. Defaults
#'        to one 1.
#' @param BPPARAM Optional `BiocParallelParam` object. If not provided, the function
#'        will configure an appropriate backend automatically.
#' @param low_cutoff Ceiling for small SD estimates. Default 0.1
#' @param high_cutoff Floor for large SD estimates. Default 10
#' @param est_disperion_prior Estimate a prior for disperion (default = FALSE)
#' @return An object of class \code{strainspy_priors}.
#' @examples
#' if (requireNamespace("fastglm", quietly = TRUE) && 
#' requireNamespace("betareg", quietly = TRUE)) {
#'   meta <- read.csv(system.file("extdata", "example_metadata.csv.gz", 
#'   package = "strainspy"))
#'   meta$Case_status <- factor(meta$Case_status)
#'   se <- read_sylph(system.file("extdata", "example_sylph_profile.tsv.gz", 
#'   package = "strainspy"), meta_data = meta)
#'   se <- filter_by_presence(se, min_nonzero = 30)
#'   pri <- compute_eb_priors(se[1:10, ], 
#'   ~ Case_status + Age_at_collection, nthreads = 1)
#'   class(pri)
#' }
#' 
#' @importFrom stats sd median
#' 
#' @export
compute_eb_priors <- function(se, design, nthreads = 1L, BPPARAM = NULL,
                              low_cutoff = 0.1, high_cutoff = 10,
                              est_disperion_prior = FALSE) {
  
  required_pkgs <- c("fastglm", "betareg")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace,
                                        logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0)
    stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))
  
  stopifnot(inherits(se, "SummarizedExperiment"))
  stopifnot(inherits(design, "formula"))
  if (!nobars_(design) == design)
    stop("Random effects not supported.")
  
  chunk_size <- 500
  
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  for (col in names(cd))
    if (is.numeric(cd[[col]]))
      cd[[col]] <- scale(cd[[col]])
  
  mx_pred <- stats::model.matrix(design, cd)
  mx_outcome <- SummarizedExperiment::assay(se)
  
  # drop unmatched
  mx_outcome <- mx_outcome[, rownames(mx_pred), drop = FALSE]
  
  # chunk rows
  n_strains <- nrow(mx_outcome)
  chunk_indices <- split(seq_len(n_strains),
                         ceiling(seq_len(n_strains) / chunk_size))
  chunk_list <- lapply(chunk_indices,
                       function(idxs) mx_outcome[idxs, , drop = FALSE])
  
  # parallel backend
  if ((nthreads > 1) & (.Platform$OS.type != "windows")) {
    if (is.null(BPPARAM))
      BPPARAM <- BiocParallel::SnowParam(workers = nthreads,
                                         progressbar = TRUE)
  } else {
    BPPARAM <- BiocParallel::SerialParam(progressbar = TRUE)
  }
  
  ## --------------------
  ## ZI priors
  ## --------------------
  
  cat("Computing fixef_zi priors...\n")
  
  rx_ZI <- BiocParallel::bplapply(chunk_list,
                                  glmBin_chunk,
                                  mx_pred = mx_pred,
                                  BPPARAM = BPPARAM)
  
  rx_ZI <- do.call(rbind, unlist(rx_ZI, recursive = FALSE))
  
  fixef_zi <- BiocParallel::bplapply(seq_len(ncol(rx_ZI)),
                                     function(i) getEst(rx_ZI, i),
                                     BPPARAM = BPPARAM)
  
  fixef_zi_med <- vapply(fixef_zi, function(x) x$med, numeric(1))
  fixef_zi_boot <- do.call(rbind,
                           lapply(fixef_zi, function(x) x$boot))
  
  ## --------------------
  ## Beta priors
  ## --------------------
  
  cat("Computing fixef priors...\n")
  
  rx_beta <- BiocParallel::bplapply(chunk_list,
                                    beta_chunk,
                                    mx_pred = mx_pred,
                                    est_disp = est_disperion_prior,
                                    BPPARAM = BPPARAM)
  
  coef_list <- unlist(lapply(rx_beta, `[[`, "coefs"), recursive = FALSE)
  
  rx_beta_mat <- do.call(rbind, coef_list)
  
  fixef <- BiocParallel::bplapply(seq_len(ncol(rx_beta_mat)),
                                  function(i) getEst(rx_beta_mat, i),
                                  BPPARAM = BPPARAM)
  
  fixef_med <- vapply(fixef, function(x) x$med, numeric(1))
  fixef_boot <- do.call(rbind,
                        lapply(fixef, function(x) x$boot))
  
  rownames(fixef_boot) <- colnames(rx_beta_mat)
  rownames(fixef_zi_boot) <- colnames(rx_ZI)
  
  ## --------------------
  ## Dispersion prior
  ## --------------------
  
  boot_disp <- matrix(numeric(0), nrow = 0, ncol = 0)
  disp_prior <- NULL
  
  if (isTRUE(est_disperion_prior)) {
    
    log_phi_vec <- unlist(lapply(rx_beta, `[[`, "log_phi"))
    log_phi_vec <- log_phi_vec[is.finite(log_phi_vec)]
    
    if (length(log_phi_vec) > 5) {
      
      B <- 1000
      boot_sd <- numeric(B)
      boot_mean <- numeric(B)
      
      for (b in seq_len(B)) {
        samp <- sample(log_phi_vec, replace = TRUE)
        boot_sd[b] <- sd(samp)
        boot_mean[b] <- mean(samp)
      }
      
      # This is likely more stable than estimating values from log_phi
      disp_sd <- median(boot_sd)
      disp_median <- median(boot_mean) 
      
      disp_prior <- data.frame(
        prior  = paste("normal(",round(disp_median, 2) , ",", round(disp_sd,2)  ,")", sep = ""),
        class    = "fixef_disp",
        coef = 1,
        stringsAsFactors = FALSE
      )
      
      boot_disp <- matrix(c(boot_mean, boot_sd), nrow = 2, byrow = TRUE)
      rownames(boot_disp) <- c("log_phi_mean", "log_phi_sd")
    } else {
      warning('Dispersion cannot be estimated, not enough strains in sample (<=5)')
    }
  }
  
  ## --------------------
  ## Build prior_df
  ## --------------------
  
  term_names <- colnames(rx_beta_mat)
  intercept_idx <- grep("(Intercept)", term_names)
  
  if (length(intercept_idx) == 1) {
    term_names <- term_names[-intercept_idx]
    sd_fixef <- fixef_med[-intercept_idx]
    sd_fixef_zi <- fixef_zi_med[-intercept_idx]
  } else {
    sd_fixef <- fixef_med
    sd_fixef_zi <- fixef_zi_med
  }
  
  prior_df <- make_prior_df(term_names,
                            sd_fixef,
                            sd_fixef_zi,
                            low_cutoff = low_cutoff,
                            high_cutoff = high_cutoff)
  
  if (!is.null(disp_prior))
    prior_df <- rbind(prior_df, disp_prior)
  
  ## --------------------
  ## Return object
  ## --------------------
  
  methods::new("strainspy_priors",
               method = "empirical",
               priors_df = prior_df,
               boot_fixef = fixef_boot,
               boot_fixef_ZI = fixef_zi_boot,
               boot_disp = boot_disp,
               design = design,
               call = match.call())
}


#' Constructor for strainspy_priors using preset or user defined priors
#' 
#' @param se A \code{SummarizedExperiment} object with count matrix and colData.
#' @param design A fixed-effects-only formula (e.g., ~ group + age + sex).
#' @param method Prior type ("preset_weak", "preset_strong", "manual"). To compute 
#' empirical priors, use `compute_eb_priors()`. Default `preset_weak`.
#' @param add_dispersion_prior Add a prior of normal(0,1) to the dispersion parameter. Default False
#' @param priors_df NULL if `method == 'preset_weak' | 'preset_strong'`, else a 
#' data.frame of priors suitable for glmmTMB model fitting. See `?glmmTMB::priors`
#' for more details. Default NULL.
#' @return An object of class \code{strainspy_priors}.
#' @examples
#' meta <- read.csv(system.file("extdata", "example_metadata.csv.gz",
#'  package = "strainspy"))
#' meta$Case_status <- factor(meta$Case_status)
#' se <- read_sylph(system.file("extdata", "example_sylph_profile.tsv.gz", 
#' package = "strainspy"), meta_data = meta)
#' pri <- suppressWarnings(define_priors(se[1:10, ], 
#' ~ Case_status + Age_at_collection, method = "preset_weak"))
#' class(pri)
#' 
#' @export
define_priors = function(se, design, method = 'preset_weak', add_dispersion_prior = FALSE, priors_df = NULL){
  stopifnot(inherits(se, "SummarizedExperiment"))
  stopifnot(inherits(design, "formula"))
  ## Check if formula has random effects
  if(!nobars_(design) == design) stop("Formula contains random effects, which are currently not supported.")
  
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  mx_pred <- stats::model.matrix(design, cd)
  
  term_names <- colnames(mx_pred)
  
  if (method %in% c("preset_weak", "preset_strong")) {
    intercept_idx = grep("(Intercept)", term_names)
    if(length(intercept_idx) == 1){
      warning("No MAP prior will be set to intercept")
      term_names = term_names[-intercept_idx]
    }
    
    priors_df <- get_preset_priors(term_names, method)
    
    if(add_dispersion_prior)
      priors_df <- rbind(priors_df, get_dispersion_prior())
  } else if(method == "manual") {
    if (is.null(priors_df)) stop("You must supply a priors_df data.frame when method = 'manual'")
    
    if(!is.data.frame(priors_df)) stop("priors_df must be a data.frame. See `?glmmTMB::priors`")
    
    required_cols <- c("prior", "class", "coef")
    missing_cols <- setdiff(required_cols, colnames(priors_df))
    if (length(missing_cols) > 0) {
      stop("Manual prior data.frame is missing required column(s): ", paste(missing_cols, collapse = ", "))
    }
    
    warning("No checks will be performed on the provided `prior`, check `?glmmTMB::priors` for details")
  } else {
    stop("Unknown method.")
  }
  
  methods::new("strainspy_priors",
               method = method,
               priors_df = priors_df,
               boot_fixef = matrix(numeric(0), nrow = 0, ncol = 0),
               boot_fixef_ZI = matrix(numeric(0), nrow = 0, ncol = 0),
               design = design,
               call = match.call())
  
}


glmBin_chunk <- function(chunk, mx_pred) {
  coef_list <- vector("list", nrow(chunk))
  keep <- logical(nrow(chunk))
  
  for (i in seq_len(nrow(chunk))) {
    y <- as.numeric(chunk[i, ] > 0)
    
    fit <- tryCatch(
      fastglm::fastglmPure(x = mx_pred, y = y, family = binomial()),
      error = function(e) NULL
    )
    
    if (!is.null(fit)) {
      coef_list[[i]] <- fit$coefficients
      keep[i] <- TRUE
    }
  }
  
  coef_list[keep]  # return only successful fits
}

beta_chunk <- function(chunk, mx_pred, est_disp = FALSE) {
  
  n <- nrow(chunk)
  
  coef_list   <- vector("list", n)
  log_phi_vec <- if (est_disp) numeric(n) else NULL
  keep        <- logical(n)
  
  for (i in seq_len(n)) {
    
    y <- chunk[i, ]
    idx <- which(y != 0)
    
    if (length(idx) < 0.1 * length(y)) next
    
    y_sub <- offset_ANI(y[idx] / 100)
    X_sub <- mx_pred[idx, , drop = FALSE]
    
    fit <- tryCatch(
      betareg::betareg.fit(X_sub, y_sub),
      error = function(e) NULL
    )
    
    if (!is.null(fit)) {
      coef_list[[i]] <- fit$coefficients$mean
      
      if (est_disp)
        log_phi_vec[i] <- fit$coefficients$precision
      
      keep[i] <- TRUE
    }
  }
  
  out <- list(
    coefs   = coef_list[keep],
    log_phi = if (est_disp) log_phi_vec[keep] else NULL
  )
  
  return(out)
}



getDispEst <- function(log_phi_vec, B = 2000, trim = 0.02, seed = NULL) {
  
  if (!is.null(seed)) {
    warning("`seed` is ignored to comply with Bioconductor policy against set.seed() in package code.")
  }
  
  x <- log_phi_vec[is.finite(log_phi_vec)]
  
  # optional trimming
  if (!is.null(trim) && trim > 0) {
    lo <- stats::quantile(x, trim)
    hi <- stats::quantile(x, 1 - trim)
    x <- x[x > lo & x < hi]
  }
  
  # point estimates
  mu_hat <- mean(x)
  sd_hat <- sd(x)
  
  # bootstrap
  boot_mu <- numeric(B)
  boot_sd <- numeric(B)
  
  for (b in seq_len(B)) {
    samp <- sample(x, replace = TRUE)
    boot_mu[b] <- mean(samp)
    boot_sd[b] <- sd(samp)
  }
  
  list(
    mu      = mu_hat,
    sd      = sd_hat,
    mu_ci   = quantile(boot_mu, c(0.025, 0.975)),
    sd_ci   = quantile(boot_sd, c(0.025, 0.975)),
    boot_mu = boot_mu,
    boot_sd = boot_sd
  )
}

getEst <- function(rx, idx){
  estV = rx[,idx] 
  boot_sd <- replicate(1000, {
    sample_est <- sample(estV, replace = TRUE)
    stats::sd(sample_est)
  })
  
  # output the median bootSD
  return(list(boot = boot_sd, med = round(stats::median(boot_sd), 2) ))
}

make_prior_df <- function(term_names, 
                          sd_fixef = rep(5, length(term_names)), 
                          sd_fixef_zi = rep(5, length(term_names)),
                          low_cutoff = 0.5,
                          high_cutoff = 10) {
  
  stopifnot(length(term_names) == length(sd_fixef), 
            length(term_names) == length(sd_fixef_zi))
  
  clamp_sd <- function(sd_vec, label) {
    original <- sd_vec
    low_idx <- which(is.na(sd_vec) | is.nan(sd_vec) | sd_vec < low_cutoff)
    high_idx <- which(sd_vec > high_cutoff)
    
    sd_vec[low_idx] <- low_cutoff
    sd_vec[high_idx] <- high_cutoff
    
    if (length(low_idx) > 0) {
      warning(sprintf(
        "Corrected %d small SD value(s) to %s",
        length(low_idx), low_cutoff
      ))
    }
    if (length(high_idx) > 0) {
      warning(sprintf(
        "Corrected %d large SD value(s) to %s",
        length(high_idx), high_cutoff
      ))
    }
    
    sd_vec
  }
  
  sd_fixef <- clamp_sd(sd_fixef, "fixef")
  sd_fixef_zi <- clamp_sd(sd_fixef_zi, "fixef_zi")
  
  data.frame(
    prior = c(
      sprintf("normal(0,%s)", sd_fixef),
      sprintf("normal(0,%s)", sd_fixef_zi)
    ),
    class = rep(c("fixef", "fixef_zi"), each = length(term_names)),
    coef  = rep(term_names, 2),
    stringsAsFactors = FALSE
  )
}

get_preset_priors = function(term_names, type = c("preset_weak", "preset_strong")) {
  type <- match.arg(type)
  sd <- if (type == "preset_weak") 5 else 1
  sd_vec <- rep(sd, length(term_names))
  
  data.frame(
    prior = c(
      sprintf("normal(0,%s)", sd_vec),
      sprintf("normal(0,%s)", sd_vec)
    ),
    class = rep(c("fixef", "fixef_zi"), each = length(term_names)),
    coef  = rep(term_names, 2),
    stringsAsFactors = FALSE
  )
}

get_dispersion_prior = function(sd){
  warning("Mean of dispersion prior is forced to 0. This is not recommended! 
          ... If you want to estimate it, use compute_eb_priors() instead")
  data.frame(
    prior = c(paste("normal(0,",sd,")", sep = "")),
    class = c("fixef_disp"),
    coef = c(1)
  )
}

