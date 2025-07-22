#' Constructor for strainspy_priors using empirical Bayes estimation
#'
#' Fits fast binomial and beta models to estimate MAP prior SDs for use in 
#' zero-inflated and ordinal beta models fitted using `GLMMTMB`. 
#' 
#' @param se A \code{SummarizedExperiment} object with count matrix and colData.
#' @param design A fixed-effects-only formula (e.g., ~ group + age + sex).
#' @param nthreads An integer specifying the number of (CPUs or workers) to use. Defaults
#'        to one 1.
#' @param BPPARAM Optional `BiocParallelParam` object. If not provided, the function
#'        will configure an appropriate backend automatically.
#' @return An object of class \code{strainspy_priors}.
#' 
#' @importFrom stats sd median
#' 
#' @export
compute_eb_priors <- function(se, design, nthreads=1L, BPPARAM=NULL) {
  # # code to read in a se with metadata
  # se <- read_sylph("../strainspy-manuscript/data/wallen/wallen_sylph_query_gtdb_220_id99.tsv.gz")
  # se <- filter_by_presence(se, min_nonzero = 72)
  # colnames(se) <- gsub("_1", "", colnames(se))
  # se = se[,which(colnames(se) != "SRR19064765" )]
  # meta <- readr::read_csv("../strainspy-manuscript/data/wallen/wallen_metadata.csv")
  # se = strainspy::modify_metadata(se, meta)
  # design = as.formula('~ Case_status + Race + Age_at_collection + Thyroid_med + Do_you_drink_alcohol')
  
  # Check for required pakcages
  required_pkgs <- c("fastglm", "betareg")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  
  if (length(missing_pkgs) > 0) {
    stop(
      "The following packages are required for empirical Bayes estimation but not installed: ",
      paste(missing_pkgs, collapse = ", "), ".\n",
      "Please install them with install.packages(c(", 
      paste0('"', missing_pkgs, '"', collapse = ", "), "))."
    )
  }
  
  stopifnot(inherits(se, "SummarizedExperiment"))
  stopifnot(inherits(design, "formula"))
  
  ## Check if formula has random effects
  if(!nobars_(design) == design) stop("Formula contains random effects, which are currently not supported.")
  
  
  chunk_size = 500 # This probably can be better tuned
  
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  
  # Scale all numeric predictors
  for (col in names(cd)) {
    if (is.numeric(cd[[col]])) cd[[col]] <- scale(cd[[col]])
  }
  
  mx_pred <- stats::model.matrix(design, cd)
  mx_outcome <- SummarizedExperiment::assay(se)
  
  # Drop unmatched samples
  dropped <- setdiff(colnames(mx_outcome), rownames(mx_pred))
  if (length(dropped) > 0) {
    warn_msg <- if (length(dropped) <= 10) {
      paste("Dropped samples:", paste(dropped, collapse = ", "))
    } else {
      sprintf("%d samples dropped. First 10: %s", length(dropped), paste(head(dropped, 10), collapse = ", "))
    }
    warning(warn_msg)
    mx_outcome <- mx_outcome[, rownames(mx_pred), drop = FALSE]
  }
  
  # Chunk the rows
  n_strains <- nrow(mx_outcome)
  chunk_indices <- split(seq_len(n_strains), ceiling(seq_len(n_strains) / chunk_size))
  chunk_list <- lapply(chunk_indices, function(idxs) mx_outcome[idxs, , drop = FALSE])
  
  # setup backend
  # Set up parallel infrastructure
  if ((nthreads > 1) & (.Platform$OS.type != "windows")) {
    # Check the operating system and set the backend accordingly
    if (is.null(BPPARAM)) {
      BPPARAM <- BiocParallel::SnowParam(workers = nthreads, progressbar = TRUE, tasks=100)
    }
  } else {
    BPPARAM <- BiocParallel::SerialParam(progressbar = TRUE)
  }
  
  # Parallel model fitting
  cat("Computing fixef_zi priors...\n")
  
  # BiocParallel::register(BiocParallel::MulticoreParam(workers = 10))  # or SnowParam on Windows
  rx_ZI <- BiocParallel::bplapply(chunk_list, glmBin_chunk, mx_pred = mx_pred, BPPARAM = BPPARAM)
  rx_ZI <- do.call(rbind, unlist(rx_ZI, recursive = FALSE))
  fixef_zi <- BiocParallel::bplapply(seq_len(ncol(rx_ZI)), function(i) getEst(rx_ZI, i), BPPARAM = BPPARAM)
  fixef_zi_med <- vapply(fixef_zi, function(x) x$med, numeric(1))
  fixef_zi_boot <- do.call(rbind, lapply(fixef_zi, function(x) x$boot))
  
  # # merge multi-levels into 1 by taking max
  # fixef_zi_idx = match_to_max_index(design, colnames(rx_ZI), fixef_zi_med)
  # fixef_zi_med = fixef_zi_med[fixef_zi_idx]
  # fixef_zi_boot = fixef_zi_boot[fixef_zi_idx, ]; rownames(fixef_zi_boot) = names(fixef_zi_idx)
  
  cat("Computing fixef priors...\n")
  rx_beta <- BiocParallel::bplapply(chunk_list, beta_chunk, mx_pred = mx_pred, BPPARAM = BPPARAM)
  rx_beta <- do.call(rbind, unlist(rx_beta, recursive = FALSE))
  fixef <- BiocParallel::bplapply(seq_len(ncol(rx_beta)), function(i) getEst(rx_beta, i), BPPARAM = BPPARAM)
  fixef_med <- vapply(fixef, function(x) x$med, numeric(1))
  fixef_boot <- do.call(rbind, lapply(fixef, function(x) x$boot))
  
  # BiocParallel::bpstop(BPPARAM)
  
  # # merge multi-levels into 1 by taking max
  # fixef_idx = match_to_max_index(design, colnames(rx_beta), fixef_med)
  # fixef_med = fixef_med[fixef_idx]
  # fixef_boot = fixef_boot[fixef_idx, ]; rownames(fixef_boot) = names(fixef_idx)
  
  # This is a last resort check that is probably unncessary
  if( all(colnames(rx_beta) == colnames(rx_ZI))){
    prior_df <- make_prior_df(colnames(rx_beta), fixef_med, fixef_zi_med)
  } else {
    stop("Mismatch is colnames between beta and bionomial fit outputs")
  }
  
  # Create object
  methods::new("strainspy_priors", 
               method = "empirical", 
               priors_df = prior_df,
               boot_fixef = fixef_boot,
               boot_fixef_ZI = fixef_zi_boot,
               design = design,
               call = match.call())
}

#' Constructor for strainspy_priors using preset or user defined priors
#' 
#' @param se A \code{SummarizedExperiment} object with count matrix and colData.
#' @param design A fixed-effects-only formula (e.g., ~ group + age + sex).
#' @param method Prior type ("preset_weak", "preset_strong", "manual"). To compute 
#' empirical priors, use `compute_eb_priors()`. Default `preset_weak`.
#' @param priors_df NULL if `method == 'preset_weak' | 'preset_strong'`, else a 
#' data.frame of priors suitable for glmmTMB model fitting. See `?glmmTMB::priors`
#' for more details. Default NULL.
#' @return An object of class \code{strainspy_priors}.
#' 
#' @export
define_priors = function(se, design, method = 'preset_weak', priors_df = NULL){
  stopifnot(inherits(se, "SummarizedExperiment"))
  stopifnot(inherits(design, "formula"))
  ## Check if formula has random effects
  if(!nobars_(design) == design) stop("Formula contains random effects, which are currently not supported.")
  
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  mx_pred <- stats::model.matrix(design, cd)
  
  term_names <- colnames(mx_pred)
  
  if (method %in% c("preset_weak", "preset_strong")) {
    priors_df <- get_preset_priors(term_names, method)
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


# multi-level effect SDs need to be merged - we'll take the SD with the max level
match_to_max_index = function(design, cnms, sds){
  t2idx = sapply(attr(terms(design), "term.labels"), function(term) {
    grep(paste0("^", term), cnms)
  }, simplify = FALSE)
  
  vapply(t2idx, function(idxs) {
    if (length(idxs) == 0) NA_integer_
    else idxs[which.max(sds[idxs])]
  }, integer(1))
  
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

beta_chunk <- function(chunk, mx_pred) {
  coef_list <- vector("list", nrow(chunk))
  keep <- logical(nrow(chunk))
  
  for (i in seq_len(nrow(chunk))) {
    y <- chunk[i, ]
    idx = which(y != 0)
    
    if (length(idx) < 0.1*nrow(mx_pred)) next
    
    y = offset_ANI(y[idx]/100)
    
    X = mx_pred[idx, ]
    
    fit <- tryCatch(
      betareg::betareg.fit(X, y),
      error = function(e) NULL
    )
    
    if (!is.null(fit)) {
      coef_list[[i]] <- fit$coefficients$mean
      keep[i] <- TRUE
    }
  }
  
  coef_list[keep]  # return only successful fits
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

