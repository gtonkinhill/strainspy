#' Simulate and Analyze Beta-Distributed Data with a Treatment Effect
#'
#' This function applies a power transformation to simulate a treatment effect.
#'
#' @param x a vector of sequence identity values to transform
#' @param lambda the power to apply to the transformation
#' @return a vector of transformed values
#' @examples
#' N_group <- 1000
#' difference <- 0.5
#'
#' shape1 <- 2
#' shape2 <- 10
#'
#' expected_difference <-  beta(shape1 + difference, shape2)/
#'      beta(shape1, shape2) - shape1/(shape1+shape2)
#'
#' df <- data.frame(cat=rep(c('a','b'), each=N_group),
#'                  vals=c(rbeta(N_group, shape1, shape2),
#'                         rescale_beta(rbeta(N_group, shape1, shape2),
#'                                      difference)))
#'
#' m <- glmmTMB::glmmTMB(vals ~ cat, data=df, family=beta_family(link="logit"))
#' coefs <- fixef(m)$cond
#'
#' # Compute predicted means on the original scale:
#' inv_logistic <- function(x) {return(exp(x)/(1+exp(x)))}
#' mu_a <- inv_logistic(coefs["(Intercept)"])
#' mu_b <- inv_logistic(coefs["(Intercept)"] + coefs["catb"])
#'
#' cat("Simulated difference on beta scale", expected_dff, "\n")
#' cat("Difference on original beta scale (group b - group a):", mu_b - mu_a, "\n")
#' @export
rescale_beta <- function(x, lambda = 0.5) {
  # Ensure x is between 0 and 1
  if (any(x < 0 | x > 1)) stop("Values must be between 0 and 1")

  # Apply power transformation (Beta regression-style scaling)
  x_rescaled <- x^lambda

  return(x_rescaled)
}





