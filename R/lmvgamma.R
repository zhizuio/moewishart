#' @title Multivariate Gamma Function
#'
#' @description
#' log of the multivariate gamma
#'
#' @name lmvgamma
#'
#' @importFrom stats rbinom rnorm runif rexp rgamma
#'
#' @param a TBA
#' @param p TBA
#'
#' @return A value of the log of the multivariate gamma
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
