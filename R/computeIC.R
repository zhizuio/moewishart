#' @title Information criteria for Wishart mixtures and MoE models
#'
#' @description
#' Compute AIC, BIC, and ICL for EM fits; and PSIS-LOO expected
#' log predictive density (elpd_loo) for Bayesian fits. Supports
#' \code{mixturewishart} (finite mixture) and \code{moewishart} (MoE with
#' covariates in gating).
#'
#' @name computeIC
#'
#' @importFrom utils tail
#' @importFrom loo loo
#'
#' @param fit A fitted object returned by \code{mixturewishart()} or
#'   \code{moewishart()}.
#'
#' @details
#' For EM fits:
#' \itemize{
#'   \item AIC: \eqn{\mathrm{AIC} = 2k - 2\ell},
#'         where \eqn{k} is the number of free parameters and
#'         \eqn{\ell} is the maximized log-likelihood (last EM iteration).
#'   \item BIC: \eqn{\mathrm{BIC} = k \log n - 2\ell}.
#'   \item ICL: \eqn{\mathrm{ICL} = \mathrm{BIC} + \sum_{i=1}^n
#'         \sum_{k=1}^K \tau_{ik}\log \tau_{ik}},
#'         i.e., BIC plus the entropy term (classification likelihood
#'         approximation).
#' }
#'
#' Parameter counting \eqn{k}:
#' \itemize{
#'   \item For \code{mixturewishart}:
#'         \eqn{k = (K-1) + K \cdot \frac{p(p+1)}{2} + K \cdot \mathbb{I}[\text{estimate }\nu]},
#'         where \eqn{(K-1)} are mixture weights, each \eqn{\Sigma_k} has
#'         \eqn{\frac{p(p+1)}{2}} free parameters, and \eqn{\nu_k} adds 1
#'         per component when estimated.
#'   \item For \code{moewishart}:
#'         \eqn{k = q\,(K-1) + K \cdot \frac{p(p+1)}{2} + K \cdot \mathbb{I}[\text{estimate }\nu]},
#'         where \eqn{q\,(K-1)} are gating regression coefficients
#'         (last class is reference with zero column).
#' }
#'
#' For Bayesian fits:
#' \itemize{
#'   \item elpd_loo: computed via \code{loo::loo} using per-observation
#'         log-likelihood draws. Requires \code{fit$loglik_individual}
#'         as a matrix of size \eqn{S \times n} (draws by observations),
#'         as produced by the provided samplers. Returns \code{elpd_loo},
#'         \code{se_elpd_loo}, \code{p_loo}, and \code{looic}.
#' }
#'
#' Notes:
#' \itemize{
#'   \item AIC/BIC/ICL are defined for MLE (EM) fits. For Bayesian fits,
#'         prefer elpd_loo (or WAIC) rather than AIC/BIC.
#'   \item The entropy term in ICL uses EM responsibilities
#'         \eqn{\tau_{ik}} (field \code{tau} for \code{mixturewishart},
#'         \code{gamma} for \code{moewishart}).
#' }
#'
#' @return
#' - For \code{method="em"}: a list with fields \code{AIC}, \code{BIC}, \code{ICL}.
#' - For \code{method="bayes"}: a list with fields ICL and elpd of class
#'   `"loo"` as returned by [loo::loo()] that contains fields `estimates`,
#'   `pointwise`, `diagnostics`
#'
#' @examples
#'
#' # Bayesian example (MoE)
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
#' set.seed(123)
#' fit <- moewishart(
#'   dat$S,
#'   X = cbind(1, dat$X), K = 3,
#'   mh_sigma = c(0.2, 0.1, 0.1), # RW-MH variances (length K)
#'   mh_beta = c(0.2, 0.2), # RW-MH variances (length K-1)
#'   niter = 100, burnin = 50
#' )
#' computeIC(fit)
#'
#' @export
computeIC <- function(fit) {
  if (inherits(fit, "mixturewishart.em") || inherits(fit, "moewishart.em")) {
    # Extract dimensions and responsibilities
    if (inherits(fit, "moewishart.em")) {
      # MoE fit returned by moewishart(..., method="em")
      n <- nrow(fit$gamma)
      K <- ncol(fit$gamma)
      p <- nrow(fit$Sigma[[1]])
      q <- nrow(fit$Beta)
      tau <- fit$gamma
      loglik <- tail(fit$loglik, 1)
      # Count parameters: gating q*(K-1), Sigma, optional nu
      k <- q * (K - 1) + K * (p * (p + 1) / 2) + fit$estimate_nu * K
    } else {
      # Mixture model fit returned by mixturewishart(..., method="em")
      n <- nrow(fit$tau)
      K <- ncol(fit$tau)
      p <- nrow(fit$Sigma[[1]])
      tau <- fit$tau
      loglik <- tail(fit$loglik, 1)
      # Count parameters: (K-1) weights, Sigma, optional nu
      k <- (K - 1) + K * (p * (p + 1) / 2) + fit$estimate_nu * K
    }
    AIC <- 2 * k - 2 * loglik
    BIC <- k * log(n) - 2 * loglik
    # ICL: BIC + sum_i sum_k tau_{ik} log tau_{ik}
    entropy_term <- sum(tau * log(tau), na.rm = TRUE)
    ICL <- BIC + entropy_term

    out <- list(AIC = AIC, BIC = BIC, ICL = ICL)
  }

  if (inherits(fit, "mixturewishart.bayes") ||
    inherits(fit, "moewishart.bayes")) {
    # Prepare MCMC output
    burnin <- fit$burnin
    niter <- fit$niter
    thin <- fit$thin
    log_lik <- fit$loglik_individual[((1 + burnin):niter) / thin, ]

    # Compute ICL
    n <- NCOL(log_lik)
    K <- dim(fit$pi_ik)[3]
    if (inherits(fit, "mixturewishart.bayes")) {
      q <- 0
    } else {
      q <- dim(fit$Beta_samples)[2]
    }
    kk <- K * (p * (p + 1) / 2 + 1) + (K - 1) * (q + 1)
    bic <- -2 * mean(log_lik) + log(n) * kk
    r_ik <- fit$pi_ik[-c(1:burnin), , ]
    r_ik <- apply(r_ik, c(2, 3), mean)
    icl <- bic - 2 * sum(r_ik * log(r_ik), na.rm = TRUE)

    out <- list()
    out$ICL <- icl

    # Use loo package
    if (!requireNamespace("loo", quietly = TRUE)) {
      stop("Package 'loo' is required. Please install.packages('loo').")
    }

    out$elpd <- loo::loo(log_lik)
  }

  return(out)
}
