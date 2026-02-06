#' @title Log multivariate gamma function
#'
#' @description
#' Compute the log of the multivariate gamma function \eqn{\log \Gamma_p(a)}
#' for dimension \eqn{p} and parameter \eqn{a}.
#'
#' @name lmvgamma
#'
#' @importFrom stats rbinom rnorm runif rexp rgamma
#'
#' @param a Numeric. Argument of \eqn{\Gamma_p(\cdot)} (often \eqn{\nu/2}
#'   in Wishart contexts).
#' @param p Integer. Dimension \eqn{p} of the multivariate gamma.
#'
#'
#' @details
#' The multivariate gamma function \eqn{\Gamma_p(a)} is defined by:
#' \deqn{
#'   \Gamma_p(a)
#'   = \pi^{\,p(p-1)/4} \prod_{j=1}^{p} \Gamma\!\left(a + \frac{1-j}{2}\right).
#' }
#'
#' Constraints: (i) \eqn{p \in \{1,2,\dots\}} (positive integer), and 
#'   (ii) \eqn{a > (p-1)/2} to keep all gamma terms finite (as in the Wishart 
#'   normalization constant).
#'
#'
#' @return A numeric scalar equal to \eqn{\log \Gamma_p(a)}.
#'
#'
#' @examples
#'
#' # simulate data
#' set.seed(123)
#' n <- 200 # subjects
#' p <- 2
#'
#' @export
lmvgamma <- function(a, p) {
  p * (p - 1) / 4 * log(pi) + sum(lgamma(a + (1 - 1:p) / 2))
}
