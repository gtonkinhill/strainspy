#' Compute the mean ANI differences for contigs with an identity effect and 
#' perform posthoc testing to establish support for each hit.
#'
#' In the zero-inflated beta model, beta effects represent changes in
#' strain identity and are modeled on the logit scale. While this scale
#' is effective for distinguishing phenotype groups statistically, the
#' resulting effect sizes are not directly interpretable as differences
#' in mean ANI. Same applies for effect sizes from the ordinal beta model.
#'
#' This function converts logit-scale beta effects into response-scale
#' mean ANI differences using marginal means and pairwise contrasts
#' (requires `emmeans`). It therefore provides an interpretable estimate of the
#' magnitude of strain-level differences between phenotype groups, along
#' with post hoc quality checks for poorly supported beta signals. 
#'
#' Posthoc testing uses two approaches. First, it counts the number of non-zero 
#' observations supporting each beta effect. Second, it compares the estimated 
#' mean ANI difference against a preset threshold. A comment is added to each 
#' hit depending on the outcome.
#'
#' This function is designed to directly work on the `tibble` returned from 
#' `top_hits()`.
#' 
#' @param se SummarizedExperiment. A `SummarizedExperiment` object containing 
#' assay data and metadata
#' @param fit A fitted `strainspy_fit`model
#' @param th `tibble` returned by `top_hits()`.
#' @param beta_min_nz Minimum non-zero proportion per phenotype group required 
#' to consider a beta hit is well supported. Hits driven by a proportion 
#' of non zero values below this threshold at least in one phenotype group are 
#' flagged `too few nonzero`. Default 0.1. Hits driven by less than or equal to 
#' 3 observations will be flagged automatically `<3 nonzero`.
#' @param beta_min_ani_diff Minimum absolute mean ANI difference required  for 
#' a beta hit to be considered meaningful. Effect sizes small than this are 
#' flagged: `small effect`. Default = 1.5e-3 corresponds to 0.15% delta ANI. 
#' @param nthreads Number of threads for parallel processing. Default `1`.
#' @param reorder_hits Logical. Return `th` after reordering hits to priorities
#' best supported ones. Default `FALSE`.
#' @param BPPARAM Optional BiocParallelParam object.
#'
#' @return A `tibble` similar to `th`, with added columns giving details on:
#' Component (Beta, ZI or both), Mean ANI_Difference, Contrast, Non-Zero Counts 
#' of contrast components, and a general `Comment` on each hit. The commetn will 
#' be `NA` if the hit passes all posthoc tests.
#' 
#' @export
#' @importFrom BiocParallel bplapply SnowParam SerialParam
comp_ani_diff_and_posthoc_test = function(se, fit, th = NULL,
                                          beta_min_nz = 0.1, 
                                          beta_min_ani_diff = 1.5e-3, 
                                          nthreads = 1L, reorder_hits = F,
                                          BPPARAM=NULL){
  
  required_pkgs <- c("emmeans")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  
  if (length(missing_pkgs) > 0) {
    stop(
      "The following packages are required for empirical Bayes estimation but not installed: ",
      paste(missing_pkgs, collapse = ", "), ".\n",
      "Please install them with install.packages(c(", 
      paste0('"', missing_pkgs, '"', collapse = ", "), "))."
    )
  }
  # # Sanity checks - th cannot be null any more
  # if (!is.null(th)) {
  #   if (!is.null(phenotype_of_interest) || !is.null(top_hit_contigs)) {
  #     stop("Provide either `th` OR (`phenotype_of_interest` and `top_hit_contigs`), not both.")
  #   }
  # } else {
  #   if (is.null(phenotype_of_interest) || is.null(top_hit_contigs)) {
  #     stop("Provide either `th` OR both `phenotype_of_interest` and `top_hit_contigs`.")
  #   }
  # }
  
  # check thresholds
  if (!is.numeric(beta_min_nz) || beta_min_nz <= 0 || beta_min_nz >= 1)
    stop("`beta_min_nz` must be between 0 and 1.")
  
  if (!is.numeric(beta_min_ani_diff) || beta_min_ani_diff < 0)
    stop("`beta_min_ani_diff` must be >= 0.")
  
  if (!inherits(se, "SummarizedExperiment"))
    stop("`se` must be a SummarizedExperiment.")
  
  if (!inherits(fit, "strainspy_fit"))
    stop("`fit` must be a strainspy_fit object.")
  
  if(is.null(fit@vcov))
    stop("`fit` must contain variance-covariance data.")
  
  # After adding support to all, reset as needed
  # if(!any(grepl("glmZiBFit", fit@call)))
  # stop("strainspy_fit must be from a call to glmZiBFit()")
  
  
  # Annotate the top hits
  th = annotate_effect(th = th, pval_thresh = attr(th, 'alpha'))
  
  if(!"beta" %in% th$component){ # no beta hits
    if(!"both" %in% th$component) {
      stop("No identity difference hits found in the provided th:
       Running `estimate_effect_sizes()` is not required.")
    }
  }
  
  top_hit_contigs = th$Contig_name[th$component !="zi"]
  
  # Could be a wrong se are the contigs present?
  missing_contigs <- setdiff(top_hit_contigs, rownames(se))
  if (length(missing_contigs) > 0) {
    stop("The following contigs are not present in `se`: ",
         paste(missing_contigs, collapse = ", "))
  } else {
    coef_idx = match(top_hit_contigs, fit@row_data$Contig_name)
  }
  
  # pick col_data from se
  col_data <- SummarizedExperiment::colData(se)
  
  # stick with scaling done during the original fit
  if (fit@scale_continuous){ 
    for (col in names(col_data)) {
      if (is.numeric(col_data[[col]])) {
        col_data[[col]] <- scale(col_data[[col]])  # Scale numeric columns
      }
    }
  }
  
  
  # Set up parallel infrastructure
  if ((nthreads > 1) & (.Platform$OS.type != "windows")) {
    # Check the operating system and set the backend accordingly
    if (is.null(BPPARAM)) {
      BPPARAM <- BiocParallel::SnowParam(workers = nthreads, progressbar = TRUE, tasks=length(top_hit_contigs))
    }
  } else {
    BPPARAM <- BiocParallel::SerialParam(progressbar = TRUE)
  }
  
  
  op_list <- BiocParallel::bplapply(
    coef_idx,
    process_contig,
    fit = fit,
    se = se,
    col_data = col_data,
    pheno = attr(th, 'phenotype_coef'),
    BPPARAM = BPPARAM
  )
  
  op <- dplyr::bind_rows(op_list)
  
  op <- op %>%
    dplyr::left_join(th, by = "Contig_name")
  
  op = posthoc_test(op, beta_min_nz, beta_min_ani_diff)
  
  if(reorder_hits){
    # Sort the hits in a more helpful way? Mixing beta + Zi is a bit of a headache
    op <- op %>%
      dplyr::arrange(
        dplyr::desc(is.na(fail_reasons)),           # valid before flagged (beta only)
        component,
        component_p
      )
  }
  
  
  return(data.frame(Contig_name = op$Contig_name,
                    Genome_file = op$Genome_file,
                    coefficient = op$coefficient,
                    std_error = op$std_error,
                    p_value = op$p_value,
                    p_adjust = op$p_adjust,
                    zi_coefficient = op$zi_coefficient,
                    zi_std_error = op$zi_std_error,
                    zi_p_value = op$zi_p_value,
                    zi_p_adjust = op$zi_p_adjust,
                    hit_component = op$component,
                    ANI_Difference = op$estimate,
                    contrast = op$contrast,
                    ref_nz_obs = op$ref_nz,
                    ref_tot_obs = op$ref_total,
                    comp_nz_obs = op$comp_nz,
                    comp_tot_obs = op$comp_total,
                    Comment = op$fail_reasons
  ))
  
}


process_contig <- function(
    i, # row_index
    fit,
    se,
    col_data,
    pheno
) {
  
  beta1  <- unlist(fit@coefficients[i,])
  # Intercept name is sometimes off (?bug)
  names(beta1) <- sub(".*[Ii]ntercept.*", "(Intercept)", names(beta1))
  lvls = c("(Intercept)", names(beta1)[pheno])
  beta1 = beta1[lvls]
  combined_formula = as.formula(paste(c("Value", as.character(fit@design)),
                                      collapse = " "))
  
  formula1 <- strainspy:::nobars_(combined_formula)
  
  V1     <- fit@vcov[[i]][lvls, lvls]
  
  cd = col_data
  cd$Value = strainspy:::offset_ANI(
    as.vector(SummarizedExperiment::assay(se)[i, ]) / 100
  )
  cd_df <- as.data.frame(cd)
  
  # subset columns to only contain variables in formula
  mf1    <- cd_df[, all.vars(combined_formula)]
  # subset to only keep relavant rows
  mdlmx = model.matrix(formula1, cd_df)
  
  ref_rows = which(mdlmx[, "(Intercept)"] == 1 & rowSums(mdlmx) == 1)
  comp_rows = which(mdlmx[, lvls[2]] == 1)
  
  # find out support
  s1 = quickly_find_support(ref_rows, cd_df)
  s2 = quickly_find_support(comp_rows, cd_df)
  
  picked_rows <- sort(c(ref_rows, comp_rows))
  
  mf1 = mf1[picked_rows, ]
  
  rg <- emmeans::qdrg(
    formula = formula1,
    data = mf1,
    coef = beta1,
    vcov = V1,
    df = Inf,
    family = fit@family
  )
  
  emm1 = emmeans::emmeans(rg, rg@model.info$terms)
  res = emmeans::contrast(emm1, method = "pairwise")
  
  data.frame(Contig_name = rownames(se)[i], ref_nz = s1[1], ref_total = s1[2], comp_nz = s2[1], comp_total = s2[2], res)
}


quickly_find_support <- function(rows, data) {
  vals <- data$Value[rows]
  c(
    nonzero = sum(vals != 0),
    # zero    = sum(vals == 0),
    total   = length(vals)
  )
}


annotate_effect <- function(top_hit_contigs = NULL, 
                            phenotype_coef = NULL, 
                            ZB_fit = NULL, 
                            pval_thresh = 0.05, 
                            th = NULL) {
  
  if (is.null(th)) {
    th <- top_hits(fit = ZB_fit, 
                   coef = phenotype_coef, 
                   alpha = 1, 
                   method = "BH") %>%
      dplyr::filter(Contig_name %in% top_hit_contigs)
  }
  
  th <- th %>%
    dplyr::mutate(
      component = dplyr::case_when(
        p_adjust < pval_thresh & zi_p_adjust < pval_thresh ~ "both",
        p_adjust < pval_thresh ~ "beta",
        zi_p_adjust < pval_thresh ~ "zi",
        TRUE ~ "none"
      ),
      component_p = dplyr::case_when(
        component == "beta" ~ p_adjust,
        component == "zi" ~ zi_p_adjust,
        component == "both" ~ pmin(p_adjust, zi_p_adjust),
        TRUE ~ NA_real_
      )
    )
  
  return(th)
}

posthoc_test <- function(op, beta_min_nz_ratio, beta_min_ani_diff){
  op <- op %>% dplyr::rowwise() %>%
    dplyr::mutate(
      fail_reasons = paste(
        c(
          if (component %in% c("beta", "both") & min(c(ref_nz, comp_nz)) <= 3) "<3 nonzero" else NULL,
          if (component %in% c("beta", "both") & min(c(ref_nz/ref_total , comp_nz/comp_total)) <= beta_min_nz_ratio) "too few nonzero" else NULL,
          if (component %in% c("beta", "both") & abs(estimate) < beta_min_ani_diff) "small effect" else NULL
        ),
        collapse = "; "
      ),
      fail_reasons = ifelse(fail_reasons == "", NA_character_, fail_reasons)
    ) %>%
    dplyr::ungroup()
  
  return(op)
}


####################
### REFIT OPTION ###
####################

#' Estimate mean ANI differences for selected contigs analysed using a Zero
#' Inflated Beta model. This function requires refitting the model and limited
#' to ZiB models. Not exported.
#' 
#' In the zero-inflated beta model, beta effects represent changes in
#' strain identity and are modeled on the logit scale. While this scale
#' is effective for distinguishing phenotype groups statistically, the
#' resulting effect sizes are not directly interpretable as differences
#' in mean ANI.
#' 
#' This function converts logit-scale beta effects into response-scale
#' mean ANI differences using marginal means and pairwise contrasts
#' (requires `emmeans`). It therefore provides an interpretable estimate of the
#' magnitude of strain-level differences between phenotype groups, along
#' with post hoc quality checks for poorly supported beta signals.
#' 
#' Contigs can be supplied directly via `top_hit_contigs`, or indirectly
#' by passing the output of `top_hits()` via `th`. If `th` is provided,
#' `phenotype_of_interest` and `top_hit_contigs` must be `NULL`.
#' 
#' @param se SummarizedExperiment. A `SummarizedExperiment` object containing
#' assay data and metadata
#' @param ZB_fit A fitted `strainspy_fit` zero-inflated beta model
#' @param th Optional `tibble` returned by `top_hits()`. If this is provided,
#' `top_hit_contigs` and `phenotype_of_interest` are extracted automatically.
#' @param phenotype_of_interest Character string specifying the phenotype.
#' Required if `th = NULL`. Default `NULL`
#' @param top_hit_contigs Character vector of contig names to analyse.
#' Required if `th = NULL`. Default `NULL`
#' @param beta_min_nz Minimum non-zero proportion per phenotype group required
#' to consider a beta hit is well supported. Hits driven by a proportion
#' of non zero values below this threshold at least in one phenotype group are
#' flagged `too few nonzero`. Default 0.1.
#' @param beta_min_ani_diff Minimum absolute mean ANI difference required  for
#' a beta hit to be considered meaningful. Effect sizes small than this are
#' flagged: `small effect`. Default is 0.01% ANI difference = 1e-4
#' @param nthreads Number of threads for parallel processing. Default `1`.
#' @param scale_continuous Logical, whether to scale numeric columns in colData.
#' Default `TRUE`
#' @param BPPARAM Optional BiocParallelParam object.
#' 
#' @return A data.frame with columns:
#'   `Contig_name`, `Contrast`, `Min_NonZero_Ratio`, `ANI_Difference`,
#'   `Significant_Component`, and `Comment`.
#' 
#' @importFrom BiocParallel bplapply SnowParam SerialParam
# estimate_effect_sizes_refit = function(se, ZB_fit, th = NULL, 
#                                        phenotype_of_interest = NULL, 
#                                        top_hit_contigs = NULL, 
#                                        beta_min_nz = 0.1, 
#                                        beta_min_ani_diff = 1e-4, 
#                                        nthreads = 1L, scale_continuous=TRUE, 
#                                        BPPARAM=NULL){
#   
#   required_pkgs <- c("emmeans")
#   missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
#   
#   if (length(missing_pkgs) > 0) {
#     stop(
#       "The following packages are required for empirical Bayes estimation but not installed: ",
#       paste(missing_pkgs, collapse = ", "), ".\n",
#       "Please install them with install.packages(c(", 
#       paste0('"', missing_pkgs, '"', collapse = ", "), "))."
#     )
#   }
#   # Sanity checks
#   if (!is.null(th)) {
#     if (!is.null(phenotype_of_interest) || !is.null(top_hit_contigs)) {
#       stop("Provide either `th` OR (`phenotype_of_interest` and `top_hit_contigs`), not both.")
#     }
#   } else {
#     if (is.null(phenotype_of_interest) || is.null(top_hit_contigs)) {
#       stop("Provide either `th` OR both `phenotype_of_interest` and `top_hit_contigs`.")
#     }
#   }
#   
#   if (!inherits(se, "SummarizedExperiment"))
#     stop("`se` must be a SummarizedExperiment.")
#   
#   if (!inherits(ZB_fit, "strainspy_fit"))
#     stop("`ZB_fit` must be a strainspy_fit object.")
#   
#   if(!any(grepl("glmZiBFit", ZB_fit@call)))
#     stop("strainspy_fit must be from a call to glmZiBFit()")
#   
#   # pick col_data from se
#   col_data <- SummarizedExperiment::colData(se)
#   
#   # pick design from fit
#   design = ZB_fit@design
#   design_terms <- all.vars(as.formula(design)) # The formula has categorical variables together
#   
#   # Repopulate variables
#   if(!is.null(th)){
#     # If th is given, phenotype is the coefficient name in ZB_fit
#     phenotype_coef = attr(th, "phenotype_coef")
#     mm <- model.matrix(ZB_fit@design, data = col_data)
#     phenotype_of_interest = design_terms[attr(mm, "assign")[phenotype_coef]]
#     top_hit_contigs = th$Contig_name
#     alpha = attr(th, "alpha")
#   } else {
#     if (!phenotype_of_interest %in% design_terms) {
#       stop(paste("Provided phenotype_of_interest:", phenotype_of_interest, "is not found in ZB_fit@design"))
#     }
#     phenotype_coef = which(design_terms == phenotype_of_interest) + 1 # Adding one for intercept
#     # is the phenotype present?
#     if (!phenotype_of_interest %in% colnames(SummarizedExperiment::colData(se))) {
#       stop("`phenotype_of_interest` not found in colData(se).")
#     }
#     
#     # This one is a bit obvious but anyway...
#     n_groups <- length(levels(SummarizedExperiment::colData(se)[[phenotype_of_interest]]))
#     
#     if (n_groups < 2) {
#       stop("`phenotype_of_interest` must be categorical with >2 groups")
#     }
#   }
#   
#   # Could be a wrong se
#   # are the contigs present?
#   missing_contigs <- setdiff(top_hit_contigs, rownames(se))
#   if (length(missing_contigs) > 0) {
#     stop("The following contigs are not present in `se`: ",
#          paste(missing_contigs, collapse = ", "))
#   }
#   
#   # check thresholds
#   if (!is.numeric(beta_min_nz) || beta_min_nz <= 0 || beta_min_nz >= 1)
#     stop("`beta_min_nz` must be between 0 and 1.")
#   
#   if (!is.numeric(beta_min_ani_diff) || beta_min_ani_diff < 0)
#     stop("`beta_min_ani_diff` must be >= 0.")
#   
#   
#   
#   # Call top hits or Annotate
#   # function(top_hit_contigs = NULL, phenotype_coef = NULL, ZB_fit = NULL, pval_thresh = 0.05, th = NULL)
#   if(is.null(th)) {
#     th = annotate_effect(top_hit_contigs = top_hit_contigs, 
#                          phenotype_coef = phenotype_coef, 
#                          ZB_fit = ZB_fit, pval_thresh = 0.05, th = NULL)
#   } else {
#     th = annotate_effect(th = th, pval_thresh = alpha)
#   }
#   
#   
#   if (scale_continuous==TRUE){
#     for (col in names(col_data)) {
#       if (is.numeric(col_data[[col]])) {
#         col_data[[col]] <- scale(col_data[[col]])  # Scale numeric columns
#       }
#     }
#   }
#   
#   
#   frml = as.formula(paste0("~", phenotype_of_interest))
#   
#   combined_formula <- as.formula(paste(c("Value", as.character(design)),
#                                        collapse = " "))
#   # We need to refit the model using only the given contigs - no significance testing - should be okay?
#   fixed_priors = ZB_fit@priors@priors_df
#   
#   # Set up parallel infrastructure
#   if ((nthreads > 1) & (.Platform$OS.type != "windows")) {
#     # Check the operating system and set the backend accordingly
#     if (is.null(BPPARAM)) {
#       BPPARAM <- BiocParallel::SnowParam(workers = nthreads, progressbar = TRUE, tasks=length(top_hit_contigs))
#     }
#   } else {
#     BPPARAM <- BiocParallel::SerialParam(progressbar = TRUE)
#   }
#   
#   op_list <- BiocParallel::bplapply(
#     top_hit_contigs,
#     process_contig_refit,
#     se = se,
#     col_data = col_data,
#     combined_formula = combined_formula,
#     design = design,
#     fixed_priors = fixed_priors,
#     frml = frml,
#     phenotype_of_interest = phenotype_of_interest,
#     BPPARAM = BPPARAM
#   )
#   
#   op <- dplyr::bind_rows(op_list)
#   
#   op <- op %>%
#     dplyr::left_join(th, by = "Contig_name")
#   
#   op = posthoc_test(op, beta_min_nz, beta_min_ani_diff)
#   
#   # Sort the hits in a more helpful way? Mixing beta + Zi is a bit of a headache
#   op <- op %>%
#     dplyr::arrange(
#       dplyr::desc(is.na(fail_reasons)),           # valid before flagged (beta only)
#       component,
#       component_p
#     )
#   
#   return(data.frame(Contig_name = op$Contig_name,
#                     Contrast = op$contrast,
#                     Min_NonZero_Ratio = op$min_nonzero_ratio,
#                     ANI_difference = op$estimate,
#                     Significant_Component = op$component,
#                     Comment = op$fail_reasons))
# }

# 
# process_contig_refit <- function(
#     contig,
#     se,
#     col_data,
#     combined_formula,
#     design,
#     fixed_priors,
#     frml,
#     phenotype_of_interest
# ) {
#   # get row index for this contig
#   row_index <- match(contig, rownames(se))
#   
#   # copy col_data to avoid side effects
#   cd <- col_data
#   cd$Value = strainspy:::offset_ANI(
#     as.vector(SummarizedExperiment::assay(se)[row_index, ]) / 100
#   )
#   
#   cd_df <- as.data.frame(cd)
#   
#   # filter out weird beta hits driven by only a few obs
#   support <- cd_df |>
#     dplyr::group_by(.data[[phenotype_of_interest]]) |>
#     dplyr::summarise(
#       nonzero_count = sum(Value > 0),
#       nonzero = sum(Value > 0)/dplyr::n(),
#       .groups = "drop"
#     )
#   
#   # fit model
#   mdl <- tryCatch(glmmTMB::glmmTMB(
#     formula = combined_formula,
#     ziformula = design,
#     data = cd_df,
#     priors = fixed_priors,
#     family = glmmTMB::beta_family(link = "logit")
#   ),
#   error = function(e) NULL)
#   if (is.null(mdl)) {
#     return(data.frame(
#       Contig_name = contig,
#       min_nonzero_ratio = min(support$nonzero),
#       min_nonzero_count = min(support$nonzero_count),
#       convergence = FALSE
#     ))
#   }
#   
#   # emmeans
#   emm <- emmeans::emmeans(mdl, frml, component = "cond") %>%
#     emmeans::regrid(transform = "response")
#   
#   # pairwise contrasts
#   res <- emmeans::contrast(emm, method = "pairwise")
#   
#   data.frame(Contig_name = contig, min_nonzero_ratio = min(support$nonzero), min_nonzero_count = min(support$nonzero_count), convergence = mdl$fit$convergence == 0 & isTRUE(mdl$sdr$pdHess), res)
# }
# 
# posthoc_test_refit <- function(op, beta_min_nz_ratio, beta_min_ani_diff){
#   op <- op %>% dplyr::rowwise() %>%
#     dplyr::mutate(
#       fail_reasons = paste(
#         c(
#           if (component %in% c("beta", "both") & !convergence) "failed convergence" else NULL,
#           if (component %in% c("beta", "both") & min_nonzero_ratio <= beta_min_nz_ratio) "too few nonzero" else NULL,
#           if (component %in% c("beta", "both") & abs(estimate) < beta_min_ani_diff) "small effect" else NULL
#         ),
#         collapse = "; "
#       ),
#       fail_reasons = ifelse(fail_reasons == "", NA_character_, fail_reasons)
#     ) %>%
#     dplyr::ungroup()
#   
#   return(op)
# }
