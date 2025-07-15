#' Constructor for strainspy_priors using empirical Bayes estimation
#'
#' Fits fast binomial and beta models to estimate MAP prior SDs for use in 
#' zero-inflated and ordinal beta models fitted using `GLMMTMB`. 
#' 
#' @param sy A \code{SummarizedExperiment} object with count matrix and colData.
#' @param design A fixed-effects-only formula (e.g., ~ group + age + sex).
#' @param chunk_size Integer. Number of strains to process per chunk. Default: 500.
#' @return An object of class \code{strainspy_priors}.
#' 
#' @importFrom stats sd median
#' 
#' @export
ebayesian_priors <- function(sy, design, chunk_size = 500) {
  
  # Check for required pakcages
  required_pkgs <- c("fastglm", "betareg", "future", "future.apply")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  
  if (length(missing_pkgs) > 0) {
    stop(
      "The following packages are required for empirical Bayes estimation but not installed: ",
      paste(missing_pkgs, collapse = ", "), ".\n",
      "Please install them with install.packages(c(", 
      paste0('"', missing_pkgs, '"', collapse = ", "), "))."
    )
  }
  
  stopifnot(inherits(sy, "SummarizedExperiment"))
  stopifnot(inherits(design, "formula"))
  
  cd <- as.data.frame(SummarizedExperiment::colData(sy))
  
  # Scale all numeric predictors
  for (col in names(cd)) {
    if (is.numeric(cd[[col]])) cd[[col]] <- scale(cd[[col]])
  }
  
  mx_pred <- stats::model.matrix(design, cd)
  mx_outcome <- SummarizedExperiment::assay(sy)
  
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
  
  # Parallel model fitting
  future::plan(future::multisession())
  
  rx_ZI <- future.apply::future_lapply(chunk_list, glmBin_chunk, mx_pred = mx_pred)
  rx_ZI <- do.call(rbind, unlist(rx_ZI, recursive = FALSE))
  fixef_zi <- unlist(future.apply::future_lapply(seq_len(ncol(rx_ZI)), function(i) round(getEst(rx_ZI, i), 2)))
  
  rx_beta <- future.apply::future_lapply(chunk_list, beta_chunk, mx_pred = mx_pred)
  rx_beta <- do.call(rbind, unlist(rx_beta, recursive = FALSE))
  fixef <- unlist(future.apply::future_lapply(seq_len(ncol(rx_beta)), function(i) round(getEst(rx_beta, i), 2)))
  
  prior_df <- make_prior_df(attr(terms(design), "term.labels"), fixef, fixef_zi)
  
  # Create object
  methods::new("strainspy_priors", 
               method = "empirical", 
               priors_df = prior_df)
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
      coef_list[[i]] <- fit$coefficients[-1]  # drop intercept
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
      coef_list[[i]] <- fit$coefficients$mean[-1]  # drop intercept
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
  return(stats::median(boot_sd))
  
  ### visualise the boot - we need this to be a separate callable function
  # hist(boot_sd, breaks = 50, col = "grey80",
  #      main = paste("Bootstrap SDs for", colnames(rx)[idx]), xlab = "SD of Estimates", border = "white")
  # 
  # # Lines for summary stats
  # abline(v = mean(boot_sd), col = "blue", lty = 2, lwd = 2)
  # abline(v = median(boot_sd), col = "red", lty = 2, lwd = 2)
  # 
  # # CI
  # ci <- quantile(boot_sd, c(0.025, 0.975))
  # abline(v = ci[1], col = "darkgreen", lty = 3, lwd = 2)
  # abline(v = ci[2], col = "darkgreen", lty = 3, lwd = 2)
  # 
  # # Legend with values
  # legend("topright", inset = 0.01,
  #        legend = c(
  #          paste0("Mean: ", round(mean(boot_sd), 4)),
  #          paste0("Median: ", round(median(boot_sd), 4)),
  #          paste0("95% CI: [", round(ci[1], 4), ", ", round(ci[2], 4), "]")
  #        ),
  #        col = c("blue", "red", "darkgreen"),
  #        lty = c(2, 2, 3), lwd = 2, bg = "white")
  
  
}

make_prior_df <- function(term_names, 
                          sd_fixef = rep(5, length(term_names)), 
                          sd_fixef_zi = rep(5, length(term_names))) {
  
  stopifnot(length(term_names) == length(sd_fixef), 
            length(term_names) == length(sd_fixef_zi))
  
  clamp_sd <- function(sd_vec, label) {
    original <- sd_vec
    low_idx <- which(is.na(sd_vec) | is.nan(sd_vec) | sd_vec < 0.01)
    high_idx <- which(sd_vec > 10)
    
    sd_vec[low_idx] <- 0.01
    sd_vec[high_idx] <- 10
    
    if (length(low_idx) > 0) {
      warning(sprintf(
        "Corrected %d %s SD value(s) < 0.01 or invalid (NA/NaN) to 0.01 at indices: %s",
        length(low_idx), label, paste(low_idx, collapse = ", ")
      ))
    }
    if (length(high_idx) > 0) {
      warning(sprintf(
        "Capped %d %s SD value(s) > 10 to 10 at indices: %s",
        length(high_idx), label, paste(high_idx, collapse = ", ")
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

# 
# 
# # Run eBayes
# eb_MAP_priors <- function(sy){
#   cd = as.data.frame(SummarizedExperiment::colData(sy)) # full meta data
#   
#   # scale numeric coloumns in mx_pred
#   for (col in names(cd)) {
#     if (is.numeric(cd[[col]])) {
#       cd[[col]] <- scale(cd[[col]])  # Scale numeric columns
#     } 
#   }
#   
#   mx_pred = stats::model.matrix(design, cd)
#   mx_outcome = SummarizedExperiment::assay(sy)
#   
#   # Get matching sample names
#   dropped <- setdiff(colnames(mx_outcome), rownames(mx_pred))
#   
#   if (length(dropped) > 0) {
#     if (length(dropped) <= 10) {
#       warning("The following samples were dropped prior to empirical Bayes fitting stage due to missing at least one value from `design` variables: ", 
#               paste(dropped, collapse = ", "))
#     } else {
#       warning(sprintf(
#         "%d samples were dropped prior to empirical Bayes fitting stage due to missing at least one value from `design` variables. First 10 dropped samples: %s",
#         length(dropped), paste(head(dropped, 10), collapse = ", ")
#       ))
#     }
#     # Subset outcome matrix to match predictor matrix
#     mx_outcome <- mx_outcome[, rownames(mx_pred), drop = FALSE]
#   }
#   
#   
#   chunk_size <- 500
#   n_strains <- nrow(mx_outcome)
#   chunk_indices <- split(seq_len(n_strains), ceiling(seq_len(n_strains)/chunk_size))
#   
#   # Pre-split outcome matrix into chunks (if sparse, subsetting preserves sparsity)
#   chunk_list <- lapply(chunk_indices, function(idxs) mx_outcome[idxs, , drop = FALSE])
#   
#   
#   plan(multisession)  # parallel plan
#   
#   rx_ZI = future_lapply(chunk_list, glmBin_chunk, mx_pred = mx_pred)
#   rx_ZI = do.call(rbind, unlist(rx_ZI, recursive = FALSE))
#   fixef_zi <- unlist(future_lapply(seq_len(ncol(rx_ZI)), function(i) {
#     round(getEst(rx_ZI, i), digits = 2)
#   }))
#   
#   rx_beta = future_lapply(chunk_list, beta_chunk, mx_pred = mx_pred)
#   rx_beta = do.call(rbind, unlist(rx_beta, recursive = FALSE))
#   fixef <- unlist(future_lapply(seq_len(ncol(rx_beta)), function(i) {
#     round(getEst(rx_beta, i), digits = 2)
#   })) 
#   
#   
#   make_prior_df(attr(terms(design), "term.labels"), fixef, fixef_zi)
# }
# 
# 
