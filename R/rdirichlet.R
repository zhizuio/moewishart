#' @title Dirichlet random sampling
#'
#' @description
#' Generate random draws from a Dirichlet distribution with parameter
#' vector \eqn{\alpha \in \mathbb{R}_+^K}. Each draw is a length-\eqn{K}
#' probability vector on the simplex.
#'
#' @name rdirichlet
#'
#' @importFrom stats rbinom rnorm runif rexp rgamma pnorm rWishart toeplitz
#'
#' @param n Integer. Number of independent Dirichlet draws to generate.
#' @param alpha Numeric vector of positive concentration parameters
#'   \eqn{\alpha = (\alpha_1,\dots,\alpha_K)}. Its length \eqn{K} defines
#'   the dimension of the simplex.
#'
#'
#' @details
#' Definition:
#'  If \eqn{Y_k \sim \mathrm{Gamma}(\alpha_k, 1)} independently for
#'         \eqn{k=1,\dots,K} (shape \eqn{\alpha_k}, rate 1), then the
#'         normalized vector \eqn{X_k = Y_k / \sum_{j=1}^K Y_j} follows
#'         \eqn{\mathrm{Dirichlet}(\alpha)}.
#'
#' Note that
#'  \code{alpha} must be a numeric vector with strictly positive entries.
#'
#'
#' @return A numeric matrix of size \eqn{n \times K}, where each row is
#'   an independent Dirichlet draw that sums to \eqn{1}.
#'
#'
#' @examples
#'
#' set.seed(123)
#' # 3-dimensional Dirichlet with asymmetric concentration
#' alpha <- c(2, 5, 3)
#' rdirichlet(5, alpha)
#'
#' @export
rdirichlet <- function(n, alpha) {
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  x / rowSums(x)
}
