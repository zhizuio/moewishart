#' @title density of Wishart distribution
#'
#' @description
#' Compute the (log) density of a \eqn{p}-dimensional Wishart distribution
#' \eqn{W_p(\nu, \Sigma)} at an SPD matrix \eqn{S}. Returns either the
#' log-density or the density depending on \code{logarithm}.
#'
#' @name dWishart
#'
#' @importFrom stats rbinom rnorm runif rexp rgamma
#'
#' @param S Numeric \eqn{p \times p} SPD matrix at which to evaluate the
#'   density.
#' @param nu Numeric. Degrees of freedom \eqn{\nu} (must exceed \eqn{p-1}).
#' @param Sigma Numeric \eqn{p \times p} SPD scale matrix \eqn{\Sigma}.
#' @param detS_val Optional numeric. Precomputed \eqn{\log|S|} to reuse;
#'   if \code{NULL}, it is computed internally.
#' @param logarithm Logical. If \code{TRUE}, return log-density;
#'   otherwise return density.
#'
#'
#' @details
#' Let \eqn{S \sim W_p(\nu, \Sigma)} with degrees of freedom \eqn{\nu}
#' and scale matrix \eqn{\Sigma} (SPD). The density is:
#' \deqn{
#'   f(S \mid \nu, \Sigma) =
#'   \frac{|S|^{(\nu - p - 1)/2}\,
#'         \exp\{-\tfrac{1}{2}\,\mathrm{tr}(\Sigma^{-1}S)\}}
#'        {2^{\nu p/2}\,|\Sigma|^{\nu/2}\,\Gamma_p(\nu/2)},
#' }
#' where \eqn{\Gamma_p(\cdot)} is the multivariate gamma function and
#' \eqn{p} is the dimension.
#'
#' Note that
#' (i) \code{detS_val} can be supplied to avoid recomputing \eqn{\log|S|},
#'  which is useful inside EM/MCMC loops, and
#'  (ii) small diagonal jitter is added internally to \eqn{S} and \eqn{\Sigma}
#'  when computing determinants or solves for numerical stability.
#'
#' Constraints: (i) \eqn{S} and \eqn{\Sigma} must be SPD, and (ii) the Wishart
#' requires \eqn{\nu > p - 1}.
#'
#'
#' @return A numeric scalar: the log-density if \code{logarithm = TRUE},
#'   otherwise the density.
#'
#' @examples
#'
#' set.seed(123)
#' p <- 3
#' # Construct an SPD Sigma
#' A <- matrix(rnorm(p * p), p, p)
#' Sigma <- crossprod(A) + diag(p) * 0.5
#' # Draw a Wishart matrix using base stats::rWishart()
#' W <- drop(rWishart(1, df = p + 5, Sigma = Sigma))
#' # Evaluate log-density at W
#' dWishart(W, nu = p + 5, Sigma = Sigma, logarithm = TRUE)
#'
#' @export
dWishart <- function(S, nu, Sigma, detS_val = NULL, logarithm = TRUE) {
  # log density of Wishart S ~ W_p(nu, Sigma)

  p <- ncol(S)

  # Calculate log-determinants safely
  if (is.null(detS_val)) {
    detS_val <- as.numeric(determinant(S + diag(p) * 1e-6, logarithm = TRUE)$modulus)
  }
  detSig_res <- determinant(Sigma + diag(p) * 1e-6, logarithm = TRUE)
  ldSigma <- as.numeric(detSig_res$modulus)

  # Calculate trace(Sigma^{-1} S) efficiently
  # solving Sigma x = S is better than explicit inverse
  invSigma_S <- sum(diag(solve(Sigma + diag(1e-6, p), S)))

  term <- ((nu - p - 1) / 2) * detS_val - (nu / 2) * ldSigma - 0.5 * invSigma_S
  ret <- term - (nu * p / 2) * log(2) - lmvgamma(nu / 2, p)

  if (!logarithm) {
    ret <- exp(ret)
  }

  return(ret)
}
