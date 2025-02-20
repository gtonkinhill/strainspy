#' Simulate and Analyze Beta-Distributed Data with a Treatment Effect
#'
#' This function applies a simple linear transformation to simulate a treatment effect.
#'
#' @param x a vector of sequence identity values to transform
#' @param lambda the power to apply to the transformation
#' @return a vector of transformed values
#' @examples
#' N_group <- 1000
#' difference <- runif(1, 0.5, 0.9)
#'
#' shape1 <- 2
#' shape2 <- 10
#'
#' sim <- rbeta(N_group, shape1, shape2)
#'
#' df <- data.frame(cat=rep(c('a','b'), each=N_group),
#'                  vals=c(sim,
#'                         rescale_beta(sim, difference)))
#'
#' m <- glmmTMB::glmmTMB(vals ~ cat, data=df,
#'                       family=glmmTMB::beta_family(link="logit"))
#' coefs <- glmmTMB::fixef(m)$cond
#'
#' cat("Simulated log-odds", log(difference), "\n",
#'     "Inferred log-odds", coefs["catb"])
#' @export
rescale_beta <- function(x, lambda = 0.5) {
  # Ensure x is between 0 and 1
  if (any(x < 0 | x > 1)) stop("Values must be between 0 and 1")
  # Ensure lambda is between 0 and 1
  if (any(lambda < 0 | lambda > 1)) stop("lambda must be between 0 and 1")

  # Apply power transformation (Beta regression-style scaling)
  x_rescaled <- x*lambda

  # Calculate expected difference to be inferred after transformation
  expected_difference <-  beta(shape1 + difference, shape2)/
       beta(shape1, shape2) - shape1/(shape1+shape2)

  return(x_rescaled)
}





