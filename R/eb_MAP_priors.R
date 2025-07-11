# mx_pred = cbind(1, as.numeric(do.call(cbind, col_data@listData)[,-1]))
# mx_outcome = SummarizedExperiment::assay(se)
# 
# chunk_size <- 500
# n_strains <- nrow(mx_outcome)
# chunk_indices <- split(seq_len(n_strains), ceiling(seq_len(n_strains)/chunk_size))
# 
# # Pre-split outcome matrix into chunks (if sparse, subsetting preserves sparsity)
# chunk_list <- lapply(chunk_indices, function(idxs) mx_outcome[idxs, , drop = FALSE])
# 
# library(future.apply)
# library(fastglm)
# 
# process_chunk <- function(chunk, mx_pred) {
#   est <- numeric(nrow(chunk))
#   keep <- logical(nrow(chunk))
#   
#   for (i in seq_len(nrow(chunk))) {
#     y <- as.numeric(chunk[i, ]>0)
#     
#     fit <- tryCatch(
#       fastglmPure(x = mx_pred, y = y, family = binomial()),
#       error = function(e) NULL
#     )
#     
#     if (!is.null(fit)) {
#       est[i] <- fit$coefficients[2]  # effect for spiked
#       keep[i] <- TRUE
#     }
#   }
#   
#   data.frame(
#     gene_idx = seq_len(nrow(chunk))[keep],
#     effect = est[keep]
#   )
# }
# 
# plan(multisession)  # parallel plan
# 
# results_list <- future_lapply(chunk_list, process_chunk, mx_pred = mx_pred)
# 
# estV = unlist(lapply(results_list, function(df) df$effect))
# 
# boot_sd <- replicate(1000, {
#   sample_est <- sample(estV, replace = TRUE)
#   sd(sample_est)
# })
# 
# ### visualise the boot
# hist(boot_sd, breaks = 50, col = "grey80",
#      main = "Bootstrap SDs", xlab = "SD of Estimates", border = "white")
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