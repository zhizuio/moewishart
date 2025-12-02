#' @title Simulate data
#'
#' @description
#' Simulate data from mixture of experts of wishart distributions
#'
#' @name simData
#'
#' @importFrom stats rbinom rnorm runif rexp rgamma pnorm rWishart toeplitz
#'
#' @param n number of observations
#' @param p dimension of the Wishart distribution
#' @param K number of latent components
#' @param Xq number of covariates for modeling subpopulation probabilities
#' @param pis vector of probabilities of subpopulations
#' @param nus vector of degrees of freedom of Wishart distributions
#' @param Sigma list of scale matrices of the Wishart distribution
#'
#' @return An object of a list
#' \itemize{
#' \item "\code{S}" - a list of covariance matrices
#' \item "\code{z}" - a vector of probabilities
#' \item "\code{Sigma_list}" - a list of scale matrices of the Wishart distribution
#' \item "\code{X}" - a matrix or NULL of covariates
#' }
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
simData <- function(n = 200, p = 2,
                    Xq = 0, K = NA,
                    pis = c(0.4, 0.6),
                    nus = c(8, 12),
                    Sigma = NULL) {
  # number of latent components
  if (Xq == 0) {
    K <- length(pis)
  } else {
    if (is.na(K)) {
      stop("Argument 'K' must be specified!")
    }
  }

  if (length(nus) != K) {
    stop("Arguments 'pis' and 'nus' must have the same length!")
  }

  # simulate covarites
  X <- NULL
  betas <- NULL
  if (Xq > 0) {
    pis <- matrix(NA, nrow = n, ncol = K)
    betas <- matrix(runif(Xq * K, -2, 2), nrow = Xq, ncol = K)
    betas[, K] <- 0
    X <- matrix(rnorm(n * Xq), nrow = n, ncol = Xq)
    expXb <- exp(X %*% betas)
    sumExpXb <- rowSums(expXb)
    pis <- apply(expXb, 2, function(xb_k) xb_k / sumExpXb)
  } else {
    pis <- matrix(pis, nrow = n, ncol = K, byrow = TRUE)
  }

  # define the scale matrix of the Wishart distribution
  if (is.null(Sigma)) {
    if (K == 2) {
      Sigma_list <- list(diag(c(1, 2)), matrix(c(2, 0.5, 0.5, 1), 2, 2))
    } else if (K == 3) {
      Sigma_list <- list(
        matrix(c(.5, .2, 0.2, 0.7), 2, 2),
        matrix(c(2.0, .6, 0.6, 1.5), 2, 2),
        matrix(c(4, .2, .2, 3), 2, 2)
      )
    } else {
      Sigma_list <- list()
      for (k in 1:K) {
        Sigma_list[[k]] <- rWishart(1, p, toeplitz((p:1) / p))[, , 1]
      }
    }
  } else {
    if (NROW(Sigma) != p || NCOL(Sigma) != p) {
      stop("The given 'Sigma' matrix has incorrect dimension!")
    }
  }

  S_list <- vector("list", n)
  z_true <- integer(n)
  for (i in 1:n) {
    k <- sample.int(K, 1, prob = pis[i, ])
    z_true[i] <- k
    # sample S ~ Wishart(nu_k, Sigma_k) using rWishart
    W <- rWishart(1, df = nus[k], Sigma = Sigma_list[[k]])[, , 1]
    S_list[[i]] <- W
  }

  list(
    S = S_list, z = z_true,
    nu = nus, pi = pis,
    Sigma_list = Sigma_list, 
    X = X, betas = betas[, -K]
  )
}
