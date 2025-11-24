#' @title Dirichlet distribution
#'
#' @description
#' TBA
#'
#' @name rdirichlet
#'
#' @importFrom stats rbinom rnorm runif rexp rgamma pnorm rWishart toeplitz
#'
#' @param n number of observations
#' @param alpha TBA
#'
#' @return An object of ...
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
rdirichlet <- function(n, alpha) {
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  x / rowSums(x)
}
