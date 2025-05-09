#' Remove random effects from model
#'
#' @param formula Formula object
#' @return A formula without mixed components
strip_random_effects <- function(formula) {
  # Convert the formula to a character string
  formula_str <- as.character(formula)

  # Split the formula string by '+'
  rhs <- formula_str[length(formula_str)]
  terms <- strsplit(rhs, "\\+")[[1]]
  terms <- trimws(terms)

  # Filter out terms containing any random effect patterns
  fixed_terms <- terms[!grepl("\\|", terms) &
                         !grepl("random\\s*\\(", terms) &
                         !grepl("re\\s*\\(", terms)]

  # Reconstruct the formula
  if (length(fixed_terms) > 1) {
    fixed_formula_str <- paste(fixed_terms, collapse = " + ")
    fixed_formula <- as.formula(paste("~", fixed_formula_str))
  } else if (length(fixed_terms) == 1) {
    fixed_formula <- as.formula(paste("~", fixed_terms))
  } else {
    fixed_formula <- as.formula("~ 1")
  }

  return(fixed_formula)
}
