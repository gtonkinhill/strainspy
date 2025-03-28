#' Remove random effects from model
#'
#' @param formula Formula object
#' @return A formula without mixed components
strip_random_effects <- function(formula) {
  # Convert the formula to a character string
  formula_str <- as.character(formula)

  # Split the formula string by '+'
  terms <- strsplit(formula_str, "\\+")[[length(formula_str)]]

  # Filter out terms containing '|' (random effects)
  fixed_terms <- terms[!grepl("\\|", terms)]

  # Reconstruct the formula
  if (length(fixed_terms) > 1) {
    fixed_formula_str <- paste(fixed_terms, collapse = " + ")
    fixed_formula <- as.formula(paste("~", fixed_formula_str))
  } else if (length(fixed_terms) == 1) {
    fixed_formula <- as.formula(paste("~", trimws(fixed_terms)))
  } else {
    fixed_formula <- as.formula("~ 1") #if nothing is left after removing random effects, return ~1
  }

  return(fixed_formula)
}

