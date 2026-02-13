#' @title Simulate data from a Wishart mixture or mixture-of-experts model
#'
#' @description
#' Generate synthetic SPD matrices from either:
#' (i) a finite mixture of Wishart components with fixed mixing proportions, or
#' (ii) a mixture-of-experts (MoE) where mixing proportions depend on covariates
#'   via a softmax gating model.
#'
#' @name simData
#'
#' @importFrom stats rbinom rnorm runif rexp rgamma pnorm rWishart toeplitz
#'
#' @param n Integer. Number of observations to simulate.
#' @param p Integer. Dimension of the Wishart distribution (matrix size
#'   \eqn{p \times p}).
#' @param K Integer. Number of latent components. Required when
#'   \code{Xq > 0}. If \code{Xq = 0}, defaults to \code{length(pis)}.
#' @param Xq Integer. Number of covariates for the gating network
#'   (MoE case). If \code{Xq = 0}, a standard mixture (no covariates)
#'   is simulated.
#' @param betas Numeric matrix \eqn{Xq \times K} of gating coefficients
#'   used when \code{Xq > 0}. If \code{NULL}, random coefficients are
#'   generated and the last column is set to zero (reference class).
#' @param pis Numeric vector of length \eqn{K} giving fixed mixture
#'   proportions when \code{Xq = 0}. Ignored when \code{Xq > 0}.
#' @param nus Numeric vector length \eqn{K}, degrees of freedom
#'   \eqn{\nu_k} for each component (must exceed \eqn{p - 1}).
#' @param Sigma Optional list length \eqn{K} of SPD scale matrices
#'   \eqn{\Sigma_k} (each \eqn{p \times p}). If \code{NULL}, defaults are
#'   generated based on \code{K} and \code{p}.
#'
#'
#' @details
#' Models:
#' \itemize{
#'   \item Fixed mixture (no covariates, \code{Xq = 0}):
#'         \eqn{z_i \sim \mathrm{Categorical}(\pi)}, and
#'         \eqn{S_i \mid z_i=k \sim W_p(\nu_k, \Sigma_k)}.
#'   \item Mixture-of-experts (covariates, \code{Xq > 0}):
#'         Let \eqn{X_i \in \mathbb{R}^{Xq}}. The mixing weights are
#'         \eqn{\pi_{ik} = \Pr(z_i=k \mid X_i)} given by softmax regression
#'         \eqn{\pi_{ik} = \exp(X_i^\top B_k) / \sum_{j=1}^K
#'         \exp(X_i^\top B_j)}. Labels \eqn{z_i} are drawn from
#'         \eqn{\mathrm{Categorical}(\pi_i)} and
#'         \eqn{S_i \mid z_i=k \sim W_p(\nu_k, \Sigma_k)}.
#' }
#'
#' Simulation steps:
#' \enumerate{
#'   \item Construct \code{pis}:
#'         \itemize{
#'           \item If \code{Xq = 0}, replicate the provided \code{pis}
#'                 over \code{n} rows.
#'           \item If \code{Xq > 0}, generate \code{X} ~ N(0, I) and compute
#'                 softmax probabilities using \code{betas} (last column set
#'                 to zero by default identifiability).
#'         }
#'   \item If \code{Sigma} is not provided, create default \eqn{\Sigma_k}
#'         matrices (SPD) depending on \code{K} and \code{p}.
#'   \item Sample labels \eqn{z_i \sim \mathrm{Categorical}(\pi_i)}.
#'   \item Draw \eqn{S_i} from \eqn{W_p(\nu_{z_i}, \Sigma_{z_i})} via
#'         \code{rWishart}.
#' }
#'
#' Note that:
#' (i) in the MoE case, no intercept is automatically added to \code{X}.
#'         Use \code{Xq} to include desired covariates; the default
#'         \code{betas} is randomly generated with \code{betas[, K] = 0}, and
#' (ii) provided \code{Sigma} must be a list of SPD \eqn{p \times p}
#'         matrices. Provided \code{nus} must exceed \eqn{p - 1}.
#'
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item \code{S}: list of length \code{n} of simulated SPD matrices
#'         \eqn{S_i}.
#'   \item \code{z}: integer vector length \code{n} of component labels.
#'   \item \code{nu}: numeric vector length \eqn{K} of degrees of freedom.
#'   \item \code{pi}: matrix \eqn{n \times K} of mixing probabilities
#'         (rows sum to \eqn{1}).
#'   \item \code{Sigma_list}: list length \eqn{K} of the scale matrices
#'         used for simulation.
#'   \item \code{X}: matrix \eqn{n \times Xq} of covariates if
#'         \code{Xq > 0}, otherwise \code{NULL}.
#'   \item \code{betas}: the gating coefficient matrix used when
#'         \code{Xq > 0}, otherwise \code{NULL}.
#' }
#'
#'
#'
#' @examples
#'
#' # simulate data from mixture model (no covariates)
#' set.seed(123)
#' n <- 200 # subjects
#' p <- 10
#' dat <- simData(n, p,
#'   K = 3,
#'   pis = c(0.35, 0.40, 0.25),
#'   nus = c(8, 12, 3)
#' )
#' str(dat)
#'
#' @export
simData <- function(n = 200, p = 2,
                    Xq = 0, K = NA,
                    betas = NULL,
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

  # check betas matrix for MoE model
  if (Xq > 0 && (!is.null(betas))) {
    if (!is.matrix(betas)) {
      stop("Argument 'betas' must be a matrix!")
    }
    dims <- dim(betas)
    if (dims[1] != Xq || dims[2] != Xq) {
      stop("Argument 'betas' has wrong dimensions!")
    }
  }

  if (Xq > 0 && is.null(betas)) {
    betas <- matrix(runif(Xq * K, -2, 2), nrow = Xq, ncol = K)
    betas[, K] <- 0
  }

  # simulate covarites
  X <- NULL
  if (Xq > 0) {
    pis <- matrix(NA, nrow = n, ncol = K)
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
    if (NROW(Sigma[[1]]) != p || NCOL(Sigma[[1]]) != p) {
      stop("The given 'Sigma' matrix has incorrect dimension!")
    }

    if (any(nus < p)) {
      stop("The given 'nus' has too small degrees of freedom!")
    }
    Sigma_list <- Sigma
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
    X = X, betas = betas
  )
}
