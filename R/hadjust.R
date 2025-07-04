#' Performs hierarchical taxonomy aware adjustment of p-values.
#'
#' This function performs multiple testing adjustment on the results of a betaGLM object
#' integrating taxonomic information. Various correction methods are available,
#' including Harmonic Mean P-value (HMP), Bonferroni and Benjamini-Hochberg (BH).
#'
#' @param object A `strainspy_fit` object.
#' @param coef The coefficient to use for p-value adjustment. Defaults to 2.
#' @param method The method for p-value adjustment. Options include "bonferroni" and Harmonic Mean P-value (HMP).
#' @param taxonomy A taxonomy data.table object. If provided, the p-values will be adjusted at each taxonomic level.
#' @param index_range description
#' @return A tibble with original and adjusted p-values.
#'
#' @importFrom tibble add_column as_tibble
#' @importFrom dplyr group_by left_join arrange mutate across all_of n rename summarise
#' @importFrom purrr map_dfr
#' @importFrom methods slot
#' @importFrom stats setNames
#'
#' @export
hadjust <- function(object, coef=2, method = "HMP", taxonomy=NULL, index_range=FALSE) {
  # Validate input
  if (!inherits(object, "strainspy_fit")) {
    stop("Input must be a strainspy_fit object.")
  }

  # Check if taxonomy data is provided
  if (is.null(taxonomy)) {
    stop("No taxonomy data provided.")
  } else if (is.character(taxonomy)) {
    tax <- read_taxonomy(taxonomy)
  } else {
    tax <- taxonomy
  }

  # Extract p-values & coefficients
  mdl = as.character(object@call)[[1]]
  if ("caseControlFit" == mdl) {
    main_model <- "Logistic"
  } else if ("glmZiBFit" == mdl | "glmFit" == mdl) {
    main_model <- "Beta"
  } else {
    main_model <- "Unknown model" # We need to expand this as we include models
  }

  beta_res <- tibble::as_tibble(slot(object, 'row_data')) |>
    tibble::add_column(coefficient=slot(object, 'coefficients')[[coef]]) |>
    tibble::add_column(std_error=slot(object, 'std_errors')[[coef]]) |>
    tibble::add_column(p_value=slot(object, 'p_values')[[coef]]) |>
    tibble::add_column(Model=main_model, .before=1)

  if (!is.null(slot(object, 'zi_coefficients'))) {
    zi_res  <- tibble::as_tibble(slot(object, 'row_data')) |>
      tibble::add_column(coefficient=slot(object, 'zi_coefficients')[[coef]]) |>
      tibble::add_column(std_error=slot(object, 'zi_std_errors')[[coef]]) |>
      tibble::add_column(p_value=slot(object, 'zi_p_values')[[coef]]) |>
      tibble::add_column(Model="Zero-Inflated", .before=1)
  }

  # Check which column in row_data is corresponds to the taxonomic data
  tax_col <- NULL
  for (i in 1:ncol(slot(object, 'row_data'))) {
    if (all(slot(object, 'row_data')[,i] %in% tax$Genome)) {
      tax_col <- i
      break
    }
  }
  if (is.null(tax_col)) {
    stop("Taxon names do not match any column found in row data.
         Taxonomy will be ignored and p-values will only be adjusted at the strain/sequence level.")
  }

  hierarchy_df <- beta_res |>
    dplyr::left_join(tax, by = setNames("Genome", names(slot(object, 'row_data'))[tax_col])) |>
    dplyr::arrange(Domain, Phylum, Class, Order, Family, Genus, Species, index) |>
    dplyr::mutate(index=1:nrow(beta_res))

  if (!is.null(slot(object, 'zi_coefficients'))){
    hierarchy_df <- rbind(
      hierarchy_df,
      zi_res |>
        dplyr::left_join(tax, by = setNames("Genome", names(slot(object, 'row_data'))[tax_col]))  |>
        dplyr::arrange(Domain, Phylum, Class, Order, Family, Genus, Species, index) |>
        dplyr::mutate(index=1:nrow(zi_res))
      )
  }

  # Remove rows which did not converge or errored
  hierarchy_df <- hierarchy_df[!is.na(hierarchy_df$coefficient),]
  # Define the grouping columns
  grouping_cols <- c("Phylum", "Class", "Order", "Family", "Genus", "Species", colnames(hierarchy_df)[[2]])
  
  if (method == "HMP") {
    L <- nrow(beta_res)
    hierarchy_df$w <- 1/nrow(beta_res)

   
    adj <- purrr::map_dfr(grouping_cols, ~{
      hierarchy_df |>
        dplyr::group_by(Model, dplyr::across(dplyr::all_of(.x))) |>  # Use all_of to refer to the grouping column
        dplyr::summarise(
          strain_count = dplyr::n(),
          mean_coefficient = mean(coefficient),
          p_adjust = harmonicmeanp::p.hmp(p_value, w, L, multilevel = TRUE) / sum(w),
          index_min = min(index),
          index_max = max(index),
          .groups = "drop"  # Drop the grouping after summarizing
        ) |>
        dplyr::rename(Name = all_of(.x)) |>  # Use all_of to rename the grouped column
        tibble::add_column(Level = .x, .before = 1)  # Add the 'Level' column before the first column
    }) |>
      dplyr::arrange(p_adjust)
  } else {
    hierarchy_df <- hierarchy_df |>
      dplyr::group_by(Model) |>
      dplyr::mutate(padj = p.adjust(p_value, method = "bonferroni")) |>
      dplyr::ungroup()

    adj <- purrr::map_dfr(grouping_cols, ~{
      hierarchy_df |>
        dplyr::group_by(Model, dplyr::across(dplyr::all_of(.x))) |>  # Use all_of to refer to the grouping column
        dplyr::summarise(
          sequence_count = dplyr::n(),
          mean_coefficient = mean(coefficient),
          p_adjust = min(padj),
          index_min = min(index),
          index_max = max(index),
          .groups = "drop"  # Drop the grouping after summarizing
        ) |>
        dplyr::rename(Name = all_of(.x)) |>  # Use all_of to rename the grouped column
        tibble::add_column(Level = .x, .before = 1)  # Add the 'Level' column before the first column
    }) |>
      dplyr::arrange(p_adjust)
  }

  if(!index_range) {
    adj$index_min <- NULL
    adj$index_max <- NULL
  }

  return(adj)
}
