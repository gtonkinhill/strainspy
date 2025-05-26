#' #' Remove random effects from model
#' #'
#' #' @param formula Formula object
#' #' @return A formula without mixed components
#' strip_random_effects <- function(formula) {
#'   # Convert the formula to a character string
#'   formula_str <- as.character(formula)
#'   
#'   # Split the formula string by '+'
#'   rhs <- formula_str[length(formula_str)]
#'   terms <- strsplit(rhs, "\\+")[[1]]
#'   terms <- trimws(terms)
#'   
#'   # Filter out terms containing any random effect patterns
#'   fixed_terms <- terms[!grepl("\\|", terms) &
#'                          !grepl("random\\s*\\(", terms) &
#'                          !grepl("re\\s*\\(", terms)]
#'   
#'   # Reconstruct the formula
#'   if (length(fixed_terms) > 1) {
#'     fixed_formula_str <- paste(fixed_terms, collapse = " + ")
#'     fixed_formula <- as.formula(paste("~", fixed_formula_str))
#'   } else if (length(fixed_terms) == 1) {
#'     fixed_formula <- as.formula(paste("~", fixed_terms))
#'   } else {
#'     fixed_formula <- as.formula("~ 1")
#'   }
#'   
#'   return(fixed_formula)
#' }

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

