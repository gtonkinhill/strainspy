#' Remove random effects from model
#' 
#' nobars_ copied from https://github.com/bbolker/reformulas 
#' This function can work with formulas like ~cohort + timepoint + (1 + timepoint | patientID), above fails
#' We can either use it like this, or add reformulas (or lme4) as a dependency
#' @param term Formula object
#' 
#' @return A formula without mixed components
nobars_ <- function(term)
{
  if (!anyBars(term)) return(term)
  if (isBar(term)) return(NULL)
  if (isAnyArgBar(term)) return(NULL)
  if (length(term) == 2) {
    nb <- nobars_(term[[2]])
    if(is.null(nb)) return(NULL)
    term[[2]] <- nb
    return(term)
  }
  nb2 <- nobars_(term[[2]])
  nb3 <- nobars_(term[[3]])
  if (is.null(nb2)) return(nb3)
  if (is.null(nb3)) return(nb2)
  term[[2]] <- nb2
  term[[3]] <- nb3
  term
}

isBar <- function(term) {
  if(is.call(term)) {
    if((term[[1]] == as.name("|")) || (term[[1]] == as.name("||"))) {
      return(TRUE)
    }
  }
  FALSE
}

isAnyArgBar <- function(term) {
  if ((term[[1]] != as.name("~")) && (term[[1]] != as.name("("))) {
    for(i in seq_along(term)) {
      if(isBar(term[[i]])) return(TRUE)
    }
  }
  FALSE
}

anyBars <- function(term) {
  any(c('|','||') %in% all.names(term))
}

#' Automatically clean contig names for plots
#' @param contig_names Vector of contig names. These will be cleaned up. This function usually works with long sylph contigs, such as those given by `top_hits(fit_ANI)$Contig_name`
#' @param max_length Numeric. Number of characters to keep (default 50)
#' @param return_map bool. Return data.frame of original and shortened contig names for mapping
#' 
#' @return A data.frame or original and shortened names, or a vector of shortened name in the same order
clean_contig_names <- function(contig_names, max_length = 50, return_map = FALSE) {
  clean_one <- function(name, max_length) {
    if (nchar(name) <= max_length) return(name)  # Skip cleaning if already short
    
    # Step 1: Remove leading accession (e.g., "CAKRLQ010000001.1")
    name <- sub("^[^ ]+\\s+", "", name)
    
    # Step 2: Extract species name (two capitalized words), fallback to isolate/bin
    species <- sub(".*?([A-Z][a-z]+\\s[a-z]+).*", "\\1", name)
    if (grepl("^whole genome", species)) species <- NA  # failed match
    
    isolate <- sub(".*isolate\\s+([\\w.-]+).*", "\\1", name, perl=TRUE)
    bin <- sub(".*bin\\.([\\w.-]+).*", "bin.\\1", name, perl=TRUE)
    
    # Step 3: Use species if found, else isolate/bin, else truncate
    if (!is.na(species) && nchar(species) <= max_length) return(species)
    if (!is.na(isolate) && nchar(isolate) <= max_length) return(paste("isolate", isolate))
    if (!is.na(bin) && nchar(bin) <= max_length) return(bin)
    
    # Step 4: Fallback: trim long names
    substr(name, 1, max_length)
  }
  
  cleaned <- unname(sapply(contig_names, function(x) clean_one(x, max_length)))
  cleaned <- ifelse(nchar(cleaned) > max_length, substr(cleaned, 1, max_length), cleaned)
  cleaned_unique <- make.unique(cleaned, sep = "_")
  
  if(return_map){
    return(data.frame(original = contig_names, short = cleaned_unique, stringsAsFactors = FALSE))
  } else {
    return(unname(cleaned_unique))
  }
  
}


#' Generate colours for Manhattan plots
#' 
#' This function generates a vector of colours for plotting, prioritising a fixed palette of 21 distinct colours.
#' If more than 21 colours are required, additional colours are generated using \code{grDevices::rainbow()}.
#'
#' @importFrom grDevices rainbow
#' 
#' @param n Integer. The number of colours needed.
#'
#' @return A character vector of \code{n} colour hex codes.
get_colors <- function(n) {
  base_colors <- c(
    '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00',
    '#cab2d6','#6a3d9a','#b15928','#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462',
    '#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5'
  )
  if (n <= length(base_colors)) {
    return(base_colors[seq_len(n)])
  } else {
    extra_needed <- n - length(base_colors)
    extra_colors <- sample(grDevices::rainbow(extra_needed))
    return(c(base_colors, extra_colors))
  }
}

#' Offset ANI by a fixed value
#' 
#' This function reduces all non-zero ANI values by `eps`
#'
#' 
#' @param x Vector of ANI values in (0,1)
#' @param eps offset value (default 1e-2)
#'
#' @return Adjusted ANI vector
offset_ANI = function(x, eps = 1e-2){
  x[x>0] <- x[x>0] - eps
  x
}

#' Resolve and construct priors for model fitting
#'
#' Ensures `MAP_prior` is interpreted correctly and returns a valid
#' `strainspy_priors` object or NULL (for unpenalised ML fits).
#'
#' @param MAP_prior One of: NULL, a character string ("preset_weak", "preset_strong"),
#'                  a prior `data.frame`, or a `strainspy_priors` object.
#' @param se A \code{SummarizedExperiment} object.
#' @param design A fixed-effects-only formula.
#'
#' @return A `strainspy_priors` object or NULL.
resolve_priors <- function(MAP_prior, se, design) {
  if (is.null(MAP_prior)) {
    return(NULL)  # unpenalised ML fit
  }
  
  if (inherits(MAP_prior, "strainspy_priors")) {
    return(MAP_prior) # it's already a strainspy_priors object
  }
  
  # We can construct it if it matches
  if (is.character(MAP_prior)) {
    if (length(MAP_prior) != 1L || is.na(MAP_prior)) {
      stop("`MAP_prior` must be a single character string when given as text; ",
           "got length ", length(MAP_prior), ".", call. = FALSE)
    }
    if (MAP_prior %in% c("preset_weak", "preset_strong")) {
      return(define_priors(se, design, method = MAP_prior))
    } else {
      stop("Character MAP_prior must be one of: 'preset_weak', 'preset_strong'.",
           call. = FALSE)
    }
  }
  
  if (is.data.frame(MAP_prior)) {
    return(define_priors(se = se, design = design, method = "manual", priors_df = MAP_prior))
  }
  
  # Anything else is a mistake, not a request for an unpenalised fit. Returning
  # NULL here would silently drop the priors and fit by maximum likelihood, so a
  # typo or a factor/list slipping through would quietly change the model. Only
  # an explicit MAP_prior = NULL selects the unpenalised fit.
  stop("`MAP_prior` must be NULL (for an unpenalised maximum-likelihood fit), ",
       "one of 'preset_weak' / 'preset_strong', a prior `data.frame`, or a ",
       "`strainspy_priors` object; got an object of class ",
       paste(sQuote(class(MAP_prior)), collapse = "/"), ".", call. = FALSE)
}

build_tax_tree = function(tax_mhp){
  tax_levels_ <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Genome")
  tax_levels_ <- intersect(tax_levels_, colnames(tax_mhp))
  
  tax_levels <- tax_levels_[
    vapply(tax_levels_, function(col) any(!is.na(tax_mhp[[col]]) & tax_mhp[[col]] != ""), logical(1))
  ]
  # strain is called genome here
  if(! ("Genome" %in% tax_levels) ){
    if("Strain" %in% tax_levels){
      tax_levels[which(tax_levels == "Strain")] = "Genome"
    } else {
      tax_levels = c(tax_levels, "Genome")
    }
  }
  
  tax_mhp[, tax_levels] <- lapply(tax_levels, function(col) as.factor(tax_mhp[[col]]))
  fml <- as.formula(paste("~", paste(tax_levels, collapse = "/")))
  
  tree <- ape::as.phylo.formula(fml, data = tax_mhp)
  return(list(tree = tree, tax_mhp = tax_mhp))
}

colour_by_tax = function(genomes, taxonomy, tax_levels = NULL){
  # This function generates a global taxonomic color palette 
  # If provided, we assume tax_levels is ordered
  tax_mhp <- taxonomy[match(genomes, taxonomy$Genome), ]
  
  # build tree for tax_mhp and update
  btt = build_tax_tree(tax_mhp)
  tree = btt$tree
  tax_mhp = btt$tax_mhp
  
  # Pre-plot to extract tip order
  p_tmp <- ggtree::ggtree(tree)
  tip_order <- p_tmp$data |>
    dplyr::filter(isTip) |>
    dplyr::arrange(y) |>
    dplyr::pull(label)
  
  if(is.null(tax_levels)){
    return(tip_order)
  }
  
  # reorder the taxonomy
  tax_mhp <- tax_mhp[match(tip_order, tax_mhp$Genome), ]; rownames(tax_mhp) <- NULL
  
  # get the phyla in order
  order_ <- unique(tax_mhp[[tax_levels[1]]])

  tax_mhp[[tax_levels[1]]] <- factor(tax_mhp[[tax_levels[1]]], levels = order_)
  
  color_palette <- setNames(
    get_colors(length(order_)),
    order_
  )
  
  return(color_palette)
}

add_tax2tophits <- function(top_hits, taxonomy, columns = c("Species")) {
  stopifnot(is.data.frame(top_hits), is.data.frame(taxonomy))
  stopifnot(is.character(columns), length(columns) >= 1)
  
  # required column for mapping
  if (!"Genome_file" %in% names(top_hits)) {
    stop("Column 'Genome_file' not found in 'top_hits'.")
  }
  if (!"Genome" %in% names(taxonomy)) {
    stop("Column 'Genome' not found in 'taxonomy'.")
  }
  
  # check requested columns exist in taxonomy
  missing_cols <- setdiff(columns, names(taxonomy))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Column(s) %s not found in 'taxonomy'.",
      paste(missing_cols, collapse = ", ")
    ))
  }
  
  # for each requested column, map values
  for (col in columns) {
    mapped <- taxonomy[[col]][match(top_hits$Genome_file, taxonomy$Genome)]
    
    # warn if any unmatched
    if (any(is.na(mapped))) {
      bad <- unique(top_hits$Genome_file[is.na(mapped)])
      warning(sprintf(
        "Column '%s': unmatched genome(s) %d (examples: %s)",
        col, length(bad), paste(utils::head(bad, 5), collapse = ", ")
      ), call. = FALSE)
    }
    
    top_hits[[col]] <- mapped
  }
  
  return(top_hits)
}

#' Identify per-feature model fits that failed
#'
#' The fitting backends signal failure in two different ways. `glmObFit()`,
#' `glmQBFit()` and `abundanceFit()` return `NULL` for a feature that could not
#' be fitted, so the failure shows up as a zero-length entry. Both `glmZiBFit()`
#' backends instead return a fully-populated list whose coefficient matrices are
#' `NA` and whose `log_likelihood` is `NA`, so the zero-length test never fires
#' there and the `NA` log-likelihood is the only signal.
#'
#' Detecting only one of the two leaves failed fits in the results, and applying
#' the `log_likelihood` test blindly is worse: `glmQBFit()` does not report a
#' log-likelihood at all, so every successful fit would look like a failure.
#' This helper handles both protocols by treating a log-likelihood that is
#' absent from *every* result as "not reported" rather than "all failed".
#'
#' It deliberately uses `vapply()` rather than `unlist()`. `unlist()` silently
#' drops the `NULL` log-likelihood of a zero-length entry, which shortens the
#' vector and shifts every later index, so the returned positions would no
#' longer line up with `results` or with `rowData()`.
#'
#' @param results A list of per-feature fit results, one element per feature,
#'   in the same order as `rowData(se)`.
#'
#' @return An increasing integer vector of positions in `results` that should be
#'   dropped, or `integer(0)` if every fit succeeded.
failed_fit_index <- function(results) {
  n <- length(results)
  if (n == 0L) return(integer(0))

  # Protocol A: the backend returned NULL for this feature.
  failed <- vapply(results, length, integer(1)) == 0L

  # Protocol B: the backend returned a populated result carrying NA log_likelihood.
  # `reports_ll` asks whether the field is present at all, which is what
  # separates "this backend does not report a log-likelihood" (glmQBFit, where
  # no result has the field) from "this fit failed" (the field is there but NA
  # or malformed).
  reports_ll <- vapply(results, function(x) !is.null(x$log_likelihood), logical(1))
  if (any(reports_ll)) {
    ll <- rep(NA_real_, n)
    ll[reports_ll] <- vapply(results[reports_ll], function(x) {
      v <- x$log_likelihood
      if (length(v) != 1L) NA_real_ else as.numeric(v)
    }, numeric(1))
    failed <- failed | (reports_ll & is.na(ll))
  }

  sort(which(failed))
}

#' Stop if sample or feature names are not unique
#'
#' `read_sylph()`, `read_metaphlan()` and `read_sourmash()` all reduce
#' `Sample_file` to a bare basename when `clean_names = TRUE`. Two input paths
#' that differ only by directory (`runA/S1.fastq.gz` and `runB/S1.fastq.gz`)
#' therefore collapse to the same name. Duplicated names are not rejected by
#' `SummarizedExperiment`, and every downstream `match()` on sample or feature
#' names would then silently attach the wrong metadata to the wrong column, so
#' fail here instead.
#'
#' @param x Character vector of names to check.
#' @param what Character. What the names identify, used in the error message
#'   (for example `"sample"` or `"feature"`).
#' @param clean_names Logical. Whether the caller applied basename cleaning; used
#'   only to tailor the suggested remedy.
#'
#' @return `x`, invisibly, if every name is unique.
check_unique_names <- function(x, what, clean_names = FALSE) {
  dup <- unique(x[duplicated(x)])
  if (length(dup) > 0L) {
    hint <- if (isTRUE(clean_names)) {
      paste0("\nThese became identical after `clean_names = TRUE` stripped the ",
             "directory and extension. Re-run with `clean_names = FALSE` to keep ",
             "the full paths, or rename the inputs so the basenames are unique.")
    } else {
      "\nRename the inputs so that each name appears only once."
    }
    stop("Duplicated ", what, " name(s) found: ",
         paste(utils::head(dup, 10), collapse = ", "),
         if (length(dup) > 10L) paste0(" (and ", length(dup) - 10L, " more)") else "",
         ".", hint, call. = FALSE)
  }
  invisible(x)
}

#' Validate a `progress` argument
#'
#' Internal helper used by the fitting functions to check that `progress` is a
#' single, non-missing logical value.
#'
#' @param progress The value supplied to a `progress` argument.
#'
#' @return `progress`, invisibly.
check_progress <- function(progress) {
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("`progress` must be a single logical value (TRUE or FALSE).", call. = FALSE)
  }
  invisible(progress)
}
