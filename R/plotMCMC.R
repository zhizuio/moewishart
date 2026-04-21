#' @title Trace plot of parameters
#'
#' @description
#' Draw MCMC trace plot of parameters
#'
#' @name plotMCMC
#'
#' @importFrom graphics legend matplot
#'
#' @param object returned object from Bayesian MoE or mixture model
#' @param estimator specified parameter for traceplot
#' @param coeff.idx choice of one covariate to show its coefficients in clusters
#'
#'
#' @return A trace plot 
#'
#'
#' @examples
#'
#' # simulate data
#' set.seed(123)
#' n <- 500 # subjects
#' p <- 2
#' # True gating coefficients (last column zero)
#' set.seed(123)
#' Xq <- 3
#' K <- 3
#' betas <- matrix(runif(Xq * K, -2, 2), nrow = Xq, ncol = K)
#' betas[, K] <- 0
#' dat <- simData(n, p,
#'   Xq = 3, K = 3, betas = betas,
#'   pis = c(0.35, 0.40, 0.25),
#'   nus = c(8, 16, 3)
#' )
#'
#' set.seed(123)
#' fit <- moewishart(
#'   dat$S,
#'   X = cbind(1, dat$X), K = 3,
#'   mh_sigma = c(0.2, 0.1, 0.1), # RW-MH variances (length K)
#'   mh_beta = c(0.2, 0.2), # RW-MH variances (length K-1)
#'   niter = 500, burnin = 200
#' )
#' 
#' plotMCMC(fit, estimator = "nu")
#'
#' @export
plotMCMC <- function(object, estimator = "nu", coeff.idx = 2) {
  
  # check if the model fitted object was from a Bayesian method
  if (inherits(object, "moewishart.bayes") || 
      inherits(object, "mixturewishart.bayes")) {
    
    if (inherits(object, "moewishart.bayes")) {
      if (!estimator %in% c("nu", "sigma", "beta", "loglik")) {
        stop("The argument 'estimator' must be one of c('nu', 'sigma', 'beta', 'loglik')!")
      }
    } else {
      if (!estimator %in% c("nu", "sigma", "pi", "loglik")) {
        stop("The argument 'estimator' must be one of c('nu', 'sigma', 'pi', 'loglik')!")
      }
    } 
  } else {
    stop("The model fitted object must be Bayesian MoE or Bayesian mixture model.")
  }
  
  # number of clusters/components
  K <- NCOL(object$nu)
  
  # traceplot for nu
  if (estimator == "nu") {
    mat <- object$nu
    matplot(mat, type = "l", lty = 1, col = 1:K, 
            xlab = "MCMC iteration", ylab = "Value",
            main = "Trace plot")
    
    legend_labels <- as.expression(
      lapply(1:K, function(i) bquote(nu[.(i)]))
    )
    legend("topright", legend = legend_labels, col = 1:K, lty = 1)
  }
  
  # traceplot for Sigma
  if (estimator == "sigma") {
    Sigma <- sapply(1:length(object$Sigma), function(xx){
      s0 <- object$Sigma[[xx]]
      sapply(1:K, function(k) { 
        as.numeric(determinant(s0[,,k], logarithm = TRUE)$modulus)
      })
    })
    mat <- t(Sigma)
    matplot(mat, type = "l", lty = 1, col = 1:K, 
            xlab = "MCMC iteration", ylab = "Value",
            main = "Trace plot")
    
    legend_labels <- as.expression(
      lapply(1:K, function(i) bquote(log~abs(Sigma[.(i)])))
    )
    legend("topright", legend = legend_labels, col = 1:K, lty = 1)
  }
  
  # traceplot for betas in MoE model
  if (estimator == "beta") {
    betas <- object$Beta_samples[, , -K] # remove (last) reference group
    coeff.idx <- 2
    mat <- betas[, coeff.idx, ] # only keep the 2nd coeff (first covariate's effect if including intercept)
    
    matplot(mat, type = "l", lty = 1, col = 1:(K-1), 
            xlab = "MCMC iteration", ylab = "Value",
            main = paste0("Trace plot of X", coeff.idx, "'s effects in various clusters"))
    
    legend_labels <- as.expression(
      lapply(1:(K-1), function(i) bquote(beta[2~.(i)]))
    )
    legend("topright", legend = legend_labels, col = 1:(K-1), lty = 1)
  }
  
  # traceplot for pi in mixture model
  if (estimator == "pi") {
    mat <- object$pi
    matplot(mat, type = "l", lty = 1, col = 1:K, 
            xlab = "MCMC iteration", ylab = "Value",
            main = "Trace plot")
    
    legend_labels <- as.expression(
      lapply(1:K, function(i) bquote(pi[.(i)]))
    )
    legend("topright", legend = legend_labels, col = 1:K, lty = 1)
  }
  
  # traceplot for log-likelihoods
  if (estimator == "loglik") {
    plot(object$loglik, type = "l", lty = 1, 
         xlab = "MCMC iteration", ylab = "Value",
         main = paste0("Trace plot of log-likelihoods"))
  }
  
}
