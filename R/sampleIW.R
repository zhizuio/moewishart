#' @title Fast sampler for Inverse-Wishart
#'
#' @description
#' TBA
#'
#' @name sampleIW
#'
#' @importFrom stats rbinom rnorm runif rexp rgamma pnorm rWishart toeplitz
#'
#' @param nu TBA
#' @param Psi_inv TBA
#'
#' @return An object of ...
#'
#' @examples
#'
#' # simulate data
#' set.seed(123)
#' n <- 200 # subjects
#' p <- 2
#'
#' @export
sampleIW <- function(nu, Psi_inv) {
  # Psi_inv is already the inverse of the scale matrix (or solved)
  # To sample Inv-Wishart(nu, Psi): sample W ~ Wishart(nu, solve(Psi))
  # Here we assume Psi_inv = solve(Psi) is passed for speed
  W <- rWishart(1, df = nu, Sigma = Psi_inv)[, , 1]
  solve(W)
}
