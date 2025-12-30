#' @title density of Wishart distribution
#'
#' @description
#' TBA
#'
#' @name dWishart
#'
#' @importFrom stats rbinom rnorm runif rexp rgamma
#'
#' @param S TBA
#' @param nu TBA
#' @param Sigma TBA
#' @param detS_val TBA
#' @param logarithm logical value for log scale or not
#'
#' @return A value of the (log) density of the Wishart distribution
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
  invSigma_S <- sum(diag(solve(Sigma, S)))

  term <- ((nu - p - 1) / 2) * detS_val - (nu / 2) * ldSigma - 0.5 * invSigma_S
  ret <- term - (nu * p / 2) * log(2) - lmvgamma(nu / 2, p)
  
  if (!logarithm) {
    ret <- exp(ret)
  }

  return(ret)
}
