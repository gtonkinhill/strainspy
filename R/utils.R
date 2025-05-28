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
clean_contig_names <- function(contig_names, max_length = 50, return_map = F) {
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


# #' Remove random effects from model
# #'
# #' @param formula Formula object
# #' @return A formula without mixed components
# strip_random_effects <- function(formula) {
#   # Convert the formula to a character string
#   formula_str <- as.character(formula)
#   
#   # Split the formula string by '+'
#   rhs <- formula_str[length(formula_str)]
#   terms <- strsplit(rhs, "\\+")[[1]]
#   terms <- trimws(terms)
#   
#   # Filter out terms containing any random effect patterns
#   fixed_terms <- terms[!grepl("\\|", terms) &
#                          !grepl("random\\s*\\(", terms) &
#                          !grepl("re\\s*\\(", terms)]
#   
#   # Reconstruct the formula
#   if (length(fixed_terms) > 1) {
#     fixed_formula_str <- paste(fixed_terms, collapse = " + ")
#     fixed_formula <- as.formula(paste("~", fixed_formula_str))
#   } else if (length(fixed_terms) == 1) {
#     fixed_formula <- as.formula(paste("~", fixed_terms))
#   } else {
#     fixed_formula <- as.formula("~ 1")
#   }
#   
#   return(fixed_formula)
# }
