#' @title Fast sampler for the inverse-Wishart distribution
#'
#' @description
#' Draw a random sample from an inverse-Wishart distribution
#' \eqn{\mathcal{IW}_p(\nu, \Psi)} using the identity
#' \eqn{S \sim \mathcal{IW}_p(\nu, \Psi)} iff
#' \eqn{S^{-1} \sim W_p(\nu, \Psi^{-1})}. This implementation accepts
#' \eqn{\Psi^{-1}} directly for speed.
#'
#' @name sampleIW
#'
#' @importFrom stats rbinom rnorm runif rexp rgamma pnorm rWishart toeplitz
#'
#' @param nu Numeric. Degrees of freedom \eqn{\nu} of the inverse-Wishart
#'   (must exceed \eqn{p - 1}).
#' @param Psi_inv Numeric \eqn{p \times p} SPD matrix equal to
#'   \eqn{\Psi^{-1}}, the inverse of the inverse-Wishart scale matrix.
#'
#'
#' @details
#' Sampling scheme:
#' \itemize{
#'   \item Sample \eqn{W \sim W_p(\nu, \Psi^{-1})} using \code{rWishart}.
#'   \item Return \eqn{S = W^{-1}}, which has
#'         \eqn{\mathcal{IW}_p(\nu, \Psi)}.
#' }
#'
#' Parameterization:
#' \itemize{
#'   \item \eqn{\nu} is the degrees of freedom, must satisfy
#'         \eqn{\nu > p - 1}.
#'   \item \eqn{\Psi} is the SPD scale matrix of the inverse-Wishart.
#'         This function expects its inverse \eqn{\Psi^{-1}} as input
#'         for efficiency (avoids repeated matrix inversions).
#' }
#'
#' Note that:
#' (i) internally calls \code{rWishart(1, df = nu, Sigma = Psi_inv)}, and
#' (ii) returns \code{solve(W)}; if numerical issues arise, consider
#'         adding a small ridge to \eqn{\Psi^{-1}} prior to sampling.
#'
#'
#' @return A \eqn{p \times p} SPD matrix \eqn{S} distributed as
#'   \eqn{\mathcal{IW}_p(\nu, \Psi)}.
#'
#' @examples
#'
#' set.seed(123)
#' p <- 3
#' # Construct an SPD scale matrix Psi
#' A <- matrix(rnorm(p * p), p, p)
#' Psi <- crossprod(A) + diag(p) * 0.5
#' Psi_inv <- solve(Psi)
#'
#' # Draw one inverse-Wishart sample
#' S <- sampleIW(nu = p + 5, Psi_inv = Psi_inv)
#' S
#'
#' @export
sampleIW <- function(nu, Psi_inv) {
  # Psi_inv is already the inverse of the scale matrix (or solved)
  # To sample Inv-Wishart(nu, Psi): sample W ~ Wishart(nu, solve(Psi))
  # Here we assume Psi_inv = solve(Psi) is passed for speed
  W <- rWishart(1, df = nu, Sigma = Psi_inv)[, , 1]
  solve(W)
}
