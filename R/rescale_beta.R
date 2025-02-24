#' Transform ANI Data to Simulate a Treatment Effect
#'
#' This function applies a simple linear transformation to simulate a treatment effect.
#'
#' @param x a vector of sequence identity values to transform
#' @param beta simulate a decreased in identiy by factor of beta (default=0.95)
#' @param zi simulate an increased shift in the proportion of zeros by zi (default=0.1)
#' @return a vector of transformed values
#' @examples
#' N_group <- 1000
#'
#' shape1 <- 2
#' shape2 <- 10
#' base_zero <- 0.5
#'
#' sim <- rbeta(N_group, shape1, shape2)
#' sim <- sim*rbinom(length(sim), 1, base_zero)
#'
#'
#' # Simulate a treatment effect.
#'
#' beta_difference <- runif(1, 0.5, 0.9)
#' treat_zero_effect <- 0.2
#'
#' treat <- rescale_beta(sim, beta_difference, treat_zero_effect)
#'
#' df <- data.frame(cat=rep(c('a','b'), each=N_group),
#'                  vals=c(sim,treat$rescaled))
#'
#' m <- glmmTMB::glmmTMB(vals ~ cat, data=df,
#'                       family=glmmTMB::beta_family(link="logit"),
#'                       ziformula = ~cat)
#'
#' coefs <- glmmTMB::fixef(m)$cond
#' zi_coefs <- glmmTMB::fixef(m)$zi
#'
#' cat("Simulated log-odds ratio", log(beta_difference), "\n",
#'     "Inferred log-odds ratio", coefs["catb"])
#'
#' cat("Simulated log-odds ratio", treat$expected_zi, "\n",
#'     "Inferred log-odds ratio", zi_coefs["catb"])
#'
#' @export
rescale_beta <- function(x, beta = 0.95, zi=0.1) {

  # Ensure x is between 0 and 1
  if (any(x < 0 | x > 1)) stop("Values must be between 0 and 1")
  # Ensure beta is between 0 and 1
  if (any(beta < 0 | beta > 1)) stop("beta must be between 0 and 1")

  # Apply power transformation (Beta regression-style scaling)
  x_rescaled <- x*beta

  p1 <- sum(x_rescaled==0)/length(x_rescaled)
  p2 <- p1 + zi

  if (p2 > 1) {
    stop("ZI is too large. The total proportion of
                   zeros is greater than 1 after adjustment")
  }

  x_rescaled[
    sample(which(x_rescaled!=0),
           size=rbinom(1, length(x_rescaled), p2) - sum(x_rescaled==0),
           replace=FALSE)] <- 0

  return(list(
    rescaled=x_rescaled,
    expected_beta=log(beta),
    expected_zi=log(p2/(1-p2)/(p1/(1-p1)))
  ))
}





