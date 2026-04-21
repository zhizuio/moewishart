#' @title EM/Bayesian estimation for Wishart MoE model
#'
#' @description
#' Fit a mixture-of-experts model for symmetric positive-definite (SPD)
#' matrices with covariate-dependent mixing proportions (gating network).
#' Components are Wishart-distributed. Supports Bayesian sampling and
#' EM-based maximum-likelihood estimation.
#'
#' @name moewishart
#'
#' @importFrom utils combn
#' @importFrom stats optim
#'
#' @param S_list List of length \eqn{n} of SPD matrices, each \eqn{p \times p}.
#'   These are the observed responses modeled by the MoE.
#' @param X Numeric matrix \eqn{n \times q} of covariates for the gating
#'   network. Include an intercept column if desired.
#' @param K Integer. Number of mixture components (experts).
#' @param niter Integer. Total iterations. Bayesian mode: total MCMC
#'   iterations (including burn-in). EM mode: maximum EM iterations.
#' @param burnin Integer. Number of burn-in iterations (Bayesian mode).
#' @param method Character; one of \code{c("bayes", "em")}. Selects
#'   sampler or optimizer.
#' @param thin Integer. Thinning interval for saving draws (Bayesian).
#' @param nu0 Numeric. Inverse-Wishart prior df for \eqn{\Sigma_k}
#'   (Bayesian). Default: \eqn{p + 2} if \code{NULL}.
#' @param Psi0 Numeric \eqn{p \times p} SPD matrix. Inverse-Wishart prior
#'   scale for \eqn{\Sigma_k} (Bayesian). Default: \code{diag(p)} if
#'   \code{NULL}.
#' @param init_nu Optional numeric vector length \eqn{K} of initial dfs
#'   \eqn{\nu_k}. Used for initialization.
#' @param estimate_nu Logical. If \code{TRUE}, estimate \eqn{\nu_k}
#'   (MH in Bayesian; Newton/EM in EM). If \code{FALSE}, keep \eqn{\nu_k}
#'   fixed at \code{init_nu}.
#' @param nu_prior_a Numeric. Prior hyperparameter \eqn{a} for \eqn{\nu_k}
#'   (Bayesian), used when \code{estimate_nu = TRUE}.
#' @param nu_prior_b Numeric. Prior hyperparameter \eqn{b} for \eqn{\nu_k}
#'   (Bayesian), used when \code{estimate_nu = TRUE}.
#' @param mh_sigma Numeric scalar or length-\eqn{K} vector. Proposal sd
#'   for MH updates on \eqn{\log(\nu_k)} (Bayesian, when estimating
#'   \eqn{\nu}).
#' @param mh_beta Numeric scalar or length-\eqn{K-1} vector. Proposal sd
#'   for MH updates of the free \eqn{B} columns (Bayesian).
#' @param sigma_beta Numeric. Prior sd of the Gaussian prior on \eqn{B}
#'   (Bayesian).
#' @param init Optional list with fields for EM initialization, e.g.,
#'   \code{beta}, \code{Sigma}, \code{nu}. See return structure.
#' @param tol Numeric. Convergence tolerance on absolute change of
#'   log-likelihood (EM), also used internally.
#' @param ridge Numeric. Small diagonal ridge added to \eqn{\Sigma_k}
#'   updates in EM for numerical stability.
#' @param verbose Logical. If \code{TRUE}, print progress information.
#'
#'
#' @details
#' MoE-Wishart Model:
#' \itemize{
#'   \item Observation: \eqn{S_i} is a \eqn{p \times p} SPD matrix. Given
#'         allocation \eqn{z_i=k}, \eqn{S_i \mid z_i \sim W_p(\nu_k,
#'         \Sigma_k)} with df \eqn{\nu_k} and scale \eqn{\Sigma_k}.
#'   \item Gating (MoE): Let \eqn{X_i} be \eqn{q}-dimensional covariates.
#'         Mixing weights \eqn{\pi_{ik} = \Pr(z_i=k \mid X_i)} follow a
#'         softmax regression:
#'         \eqn{\pi_{ik} = \exp(\eta_{ik})/\sum_{j=1}^K \exp(\eta_{ij})},
#'         where \eqn{\eta_i = X_i^\top B}, \eqn{B} is
#'         \eqn{q \times K}. Identifiability: last column of \eqn{B}
#'         is fixed to zero.
#' }
#'
#' Algorithms:
#' \enumerate{
#'   \item Bayesian (\code{method = "bayes"}): Metropolis-within-Gibbs
#'         sampler for \eqn{z}, \eqn{\Sigma_k}, optional \eqn{\nu_k}, and
#'         \eqn{B}. Gaussian priors on \eqn{B} with sd
#'         \code{sigma_beta}. Proposals use \code{mh_sigma} for
#'         \eqn{\log(\nu_k)} and \code{mh_beta} for \eqn{B}.
#'   \item EM (\code{method = "em"}): E-step responsibilities using
#'         Wishart log-densities and softmax gating. M-step updates
#'         \eqn{\Sigma_k}, optional \eqn{\nu_k}, and \eqn{B} via
#'         weighted multinomial logistic regression (BFGS).
#' }
#'
#' Note that:
#' (i) include an intercept column in \code{X}; none is added by default, and
#' (ii) all \code{S_list} elements must be SPD. A small \code{ridge} may be
#' added for stability.
#'
#'
#' @return A list whose fields depend on \code{method}:
#' \itemize{
#'   \item For \code{method = "bayes"}:
#'     \itemize{
#'       \item \code{Beta_samples}: array (\code{nsave} x \code{q} x
#'             \code{K}), saved draws of \eqn{B} (last column zero).
#'       \item \code{nu}: matrix (\code{nsave} x \code{K}), draws
#'             of \eqn{\nu_k}.
#'       \item \code{Sigma}: list of length \code{nsave}; each
#'             element is an array (\eqn{p \times p \times K}) of
#'             \eqn{\Sigma_k} draws.
#'       \item \code{z_samples}: matrix (\code{nsave} x \code{n}), draws
#'             of allocations.
#'       \item \code{pi_ik}: array (\code{nsave} x \code{n} x \code{K}),
#'             per-observation gating probabilities.
#'       \item \code{pi_mean}: matrix (\code{n} x \code{K}), posterior
#'             mean of gating probabilities.
#'       \item \code{loglik}: numeric vector (length \code{niter}),
#'             log-likelihood trace.
#'       \item \code{loglik_individual}: matrix (\code{niter} x
#'             \code{n}), per-observation log-likelihood.
#'     }
#'   \item For \code{method = "em"}:
#'     \itemize{
#'       \item \code{K, p, q, n}: problem dimensions.
#'       \item \code{Beta}: matrix (\eqn{q \times K}), gating coefficients
#'             with last column zero (reference class).
#'       \item \code{Sigma}: list length \code{K}, each a \eqn{p \times p}
#'             SPD matrix (scale).
#'       \item \code{nu}: numeric vector length \code{K}, degrees of
#'             freedom.
#'       \item \code{gamma}: matrix (\eqn{n \times K}), final
#'             responsibilities.
#'       \item \code{loglik}: numeric vector, log-likelihood by EM
#'             iteration.
#'       \item \code{iter}: integer, number of EM iterations performed.
#'     }
#' }
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
#' # Posterior means for degrees of freedom of Wishart distributions:
#' nu_mcmc <- fit$nu[-c(1:fit$burnin), ]
#' colMeans(nu_mcmc)
#'
#' @export
moewishart <- function(S_list,
                       X, # n x q matrix of covariates for gating
                       K,
                       niter = 3000,
                       burnin = 1000,
                       method = "bayes",
                       thin = 1,
                       nu0 = NULL,
                       Psi0 = NULL,
                       init_nu = NULL,
                       estimate_nu = TRUE,
                       nu_prior_a = 2, nu_prior_b = 0.1,
                       mh_sigma = 0.1,
                       mh_beta = 0.05, # MH proposal sd for gating coeffs
                       sigma_beta = 10, # Gaussian prior sd for beta
                       init = NULL,
                       tol = 1e-6,
                       ridge = 1e-8,
                       verbose = TRUE) {
  # Mixture-of-Experts Gibbs sampler for Wishart clusters
  #   S_list: list of n SPD matrices (p x p)
  #   X      : n x q covariate matrix for gating network (include intercept if desired)
  #   K      : number of experts/clusters
  #
  # Returns samples for Sigma_k, nu_k, z, and gating coefficients Beta


  if (!method %in% c("bayes", "em")) {
    stop("Argument 'method' must be either 'bayes' or 'em'!")
  }

  # TODO: remove redundant code for the common calculations between full Bayesian and EM algorithm
  if (method == "bayes") {
    n <- length(S_list)
    p <- nrow(S_list[[1]])
    q <- ncol(X)
    if (n != nrow(X)) stop("Number of rows in X must equal length(S_list).")

    if (length(mh_sigma) != K && length(mh_sigma) != 1) {
      stop("Argument 'mh_sigma' must have length 1 or K!")
    } else if (length(mh_sigma) == 1) {
      mh_sigma <- rep(mh_sigma, K)
    }
    if (length(mh_beta) != K - 1 && length(mh_beta) != 1) {
      stop("Argument 'mh_beta' must have length 1 or K-1!")
    } else if (length(mh_beta) == 1) {
      mh_beta <- rep(mh_beta, K - 1)
    }

    # Priors / defaults
    if (is.null(nu0)) nu0 <- p + 2
    if (is.null(Psi0)) Psi0 <- diag(p)
    if (is.null(init_nu)) init_nu <- rep(p + 2, K)

    # Vectorize S for fast trace computations
    S_mat <- t(sapply(S_list, as.vector)) # n x p^2
    log_det_S <- sapply(S_list, function(x) as.numeric(determinant(x, logarithm = TRUE)$modulus))

    # initialize z (kmeans on vectorized matrices)
    km <- kmeans(S_mat, centers = K, nstart = 5)
    z <- km$cluster

    # initialize cluster params
    Sigma_k <- array(0, c(p, p, K))
    nu_k <- init_nu
    for (k in 1:K) {
      idx <- which(z == k)
      if (length(idx) == 0) {
        Sigma_k[, , k] <- Psi0 / (nu0 - p - 1 + 1e-8)
      } else {
        S_sum_vec <- colSums(S_mat[idx, , drop = FALSE])
        S_sum <- matrix(S_sum_vec, p, p)
        Sigma_k[, , k] <- (Psi0 + S_sum) / (nu0 + length(idx) * nu_k[k] - p - 1)
      }
    }

    # Initialize gating coefficients Beta: q x (K-1), last column zero for identifiability
    Beta <- matrix(0, nrow = q, ncol = K) # we keep full K but enforce Beta[,K] = 0
    Beta[, 1:(K - 1)] <- matrix(rnorm(q * (K - 1), 0, 0.1), nrow = q, ncol = K - 1)
    Beta[, K] <- 0

    pi_ik <- compute_pi_ik(X, Beta) # n x K

    # Storage
    ## nsave <- floor((niter - burnin) / thin)
    nsave <- floor(niter / thin)
    if (nsave < 1) nsave <- 1
    out_beta <- array(NA, dim = c(nsave, q, K)) # store Beta at saves
    out_nu <- matrix(NA, nrow = nsave, ncol = K)
    out_Sigma <- vector("list", nsave)
    out_z <- matrix(NA, nrow = nsave, ncol = n)
    out_pi_ik <- array(NA, dim = c(nsave, n, K))
    out_pi_mean <- array(0, dim = c(n, K)) # accumulate posterior mean of pi_ik
    logliks <- numeric(niter)
    logliks_individual <- matrix(NA, nrow = niter, ncol = n)
    iter_save <- 0

    # pre-alloc
    logpost <- matrix(0, n, K)
    acc_count_nu <- numeric(K) # record acceptance counts of MH for nu
    acc_count_beta <- c(numeric(K - 1), NA) # record acceptance counts of MH for beta
    free_idx_cols <- 1:(K - 1)

    # Spike-and-slab indicators for Beta (q x (K-1))
    # pi_gamma <- 0.5
    # Gamma <- matrix(rbinom(q*(K-1), 1, pi_gamma), nrow = q, ncol = K-1)
    # Gama_samples =  array(NA, dim = c(nsave, q, K-1))

    # start_time <- Sys.time()
    for (iter in 1:niter) {
      # --- Step 1: compute log-likelihood parts per cluster (vectorized) ---
      for (k in 1:K) {
        Sig <- Sigma_k[, , k]
        chol_Sig <- tryCatch(chol(Sig), error = function(e) chol(Sig + diag(1e-8, p)))
        log_det_Sig <- 2 * sum(log(diag(chol_Sig)))
        Sig_inv <- chol2inv(chol_Sig)
        tr_val <- S_mat %*% as.vector(Sig_inv) # n x 1
        nu <- nu_k[k]
        term1 <- (nu - p - 1) / 2 * log_det_S
        term2 <- -0.5 * tr_val
        term3 <- -(nu * p / 2) * log(2) - (nu / 2) * log_det_Sig
        term4 <- -lmvgamma(nu / 2, p)
        # log posterior (pointwise) = log pi_ik + likelihood terms
        logpost[, k] <- log(pi_ik[, k] + 1e-300) + term1 + term2 + term3 + term4
      }

      # --- Step 2: sample z (categorical per i using pi_ik and likelihood) ---
      for (i in 1:n) {
        lp <- logpost[i, ]
        lp <- lp - max(lp)
        prob <- exp(lp)
        prob <- prob / sum(prob)
        z[i] <- sample.int(K, 1, prob = prob)
      }

      # --- Step 3: update gating coefficients Beta via MH ---
      # We'll update the (K-1) free columns jointly (size q*(K-1)).
      # Flatten current free parameters

      # free_idx_cols <- 1:(K - 1)

      # Update Beta columnwise
      for (free_idx_cols in 1:(K - 1)) {
        Beta_free <- as.vector(Beta[, free_idx_cols]) # length q*(K-1)
        prop <- Beta_free + rnorm(length(Beta_free), 0, mh_beta[free_idx_cols])
        Beta_prop <- Beta
        Beta_prop[, free_idx_cols] <- prop # matrix(prop, nrow = q, ncol = K - 1)
        # Beta_prop[, K] <- 0
        # compute new pi_ik (n x K)
        pi_prop <- compute_pi_ik(X, Beta_prop)
        # log prior for Beta (Gaussian iid)
        lp_prior_old <- -0.5 * sum(Beta^2) / (sigma_beta^2)
        lp_prior_new <- -0.5 * sum(Beta_prop^2) / (sigma_beta^2)
        # log-likelihood of labels given gating: sum_i log pi_i,z[i] (only depends on z)
        ll_old <- sum(log(pi_ik[cbind(1:n, z)] + 1e-300))
        ll_new <- sum(log(pi_prop[cbind(1:n, z)] + 1e-300))
        log_accept <- (ll_new + lp_prior_new) - (ll_old + lp_prior_old)
        if (log(runif(1)) < log_accept) {
          Beta <- Beta_prop
          pi_ik <- pi_prop
          # optionally track acceptance
          acc_count_beta[free_idx_cols] <- acc_count_beta[free_idx_cols] + 1
        }
      }

      # Identifiability constraint
      Beta[, K] <- 0


      # # --- Step 3b: update Gamma indicators (Gibbs) ---
      # # use exact posterior p(gamma=1 | beta) \propto pi_gamma * N(beta | 0, sigma_beta^2)
      # beta_cur <- Beta[, free_idx_cols]  # q x (K-1)
      # for (j in 1:q) {
      #   for (k in 1:(K-1)) {
      #     bval <- beta_cur[j, k]
      #     # log-prob (unnormalized) for gamma = 1 and gamma = 0
      #     log_p1 <- log(pi_gamma + 1e-300) - 0.5 * (bval^2) / (sigma_beta^2) - 0.5 * log(2 * pi * sigma_beta^2)
      #     log_p0 <- log(1 - pi_gamma + 1e-300) - 0.5 * (bval^2) / (tau0^2) - 0.5 * log(2 * pi * tau0^2)
      #     # numeric stable conversion to probability
      #     maxlog <- max(log_p1, log_p0)
      #     p1 <- exp(log_p1 - maxlog)
      #     p0 <- exp(log_p0 - maxlog)
      #     prob1 <- p1 / (p1 + p0)
      #     Gamma[j, k] <- rbinom(1, 1, prob1)
      #   }
      # }

      # --- Step 4: update Sigma_k using current assignments z ---
      for (k in 1:K) {
        idx <- which(z == k)
        nk <- length(idx)
        if (nk == 0) {
          # sample from prior-ish fallback
          Sigma_k[, , k] <- Psi0 / (nu0 - p - 1 + 1e-8)
        } else {
          S_sum_vec <- colSums(S_mat[idx, , drop = FALSE])
          S_sum <- matrix(S_sum_vec, p, p)
          nu_post <- nu0 + nk * nu_k[k]
          Psi_post <- Psi0 + S_sum
          Psi_post_inv <- tryCatch(solve(Psi_post), error = function(e) solve(Psi_post + diag(1e-8, p)))
          Sigma_k[, , k] <- sampleIW(nu_post, Psi_post_inv)
        }
      }

      # --- Step 5: update nu_k (MH) as before ---
      if (estimate_nu) {
        for (k in 1:K) {
          curr_nu <- nu_k[k]
          prop_log <- rnorm(1, log(curr_nu), mh_sigma[k])
          prop_nu <- exp(prop_log)
          if (prop_nu > p - 1 + 1e-8) {
            idx <- which(z == k)
            if (length(idx) > 0) {
              sum_log_det_S_k <- sum(log_det_S[idx])
              Sig <- Sigma_k[, , k]
              chol_Sig <- tryCatch(chol(Sig), error = function(e) chol(Sig + diag(1e-8, p)))
              log_det_Sig <- 2 * sum(log(diag(chol_Sig)))
              # Sig_inv <- chol2inv(chol_Sig)
              # S_sum_vec <- colSums(S_mat[idx, , drop = FALSE])
              # sum_tr_val <- sum(as.vector(Sig_inv) * S_sum_vec)
              calc_ll_nu <- function(val_nu) {
                term1 <- (val_nu - p - 1) / 2 * sum_log_det_S_k
                # term2 <- -0.5 * sum_tr_val ## this term is indep. of nu_k
                term3 <- -length(idx) * ((val_nu * p / 2) * log(2) + (val_nu / 2) * log_det_Sig)
                term4 <- -length(idx) * lmvgamma(val_nu / 2, p)
                term1 + term3 + term4
              }
              ll_old <- calc_ll_nu(curr_nu)
              ll_new <- calc_ll_nu(prop_nu)
            } else {
              ll_old <- 0
              ll_new <- 0
            }
            lp_old <- (nu_prior_a - 1) * log(curr_nu) - nu_prior_b * curr_nu + log(curr_nu)
            lp_new <- (nu_prior_a - 1) * log(prop_nu) - nu_prior_b * prop_nu + log(prop_nu)
            if (log(runif(1)) < (ll_new + lp_new) - (ll_old + lp_old)) {
              nu_k[k] <- prop_nu
              acc_count_nu[k] <- acc_count_nu[k] + 1
            }
          }
        }
      }

      # --- Compute loglik approx for monitoring ---
      max_l <- apply(logpost, 1, max)
      row_sums <- exp(logpost - max_l)
      logliks[iter] <- sum(max_l + log(rowSums(row_sums)))
      logliks_individual[iter, ] <- max_l + log(rowSums(row_sums))

      # --- Save samples after burnin and thinning ---
      ## if (iter > burnin && ((iter - burnin) %% thin == 0)) {
      if (iter %% thin == 0) {
        iter_save <- iter_save + 1
        if (iter_save <= nsave) {
          out_beta[iter_save, , ] <- Beta
          out_nu[iter_save, ] <- nu_k
          out_Sigma[[iter_save]] <- Sigma_k
          out_z[iter_save, ] <- z
          out_pi_ik[iter_save, , ] <- pi_ik
        }
        out_pi_mean <- out_pi_mean + pi_ik
      }

      if (verbose && (iter %% 500 == 0 || iter == 1)) {
        # elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        # rate <- iter / elapsed
        cat(sprintf(
          "Iter %4d | LL=%.1f | acc_rate_nu=%.3f | acc_rate_beta=%.3f\n",
          iter, logliks[iter], acc_count_nu / iter, # min(acc_count_nu / iter), max(acc_count_nu / iter),
          acc_count_beta / iter # min(acc_count_beta / iter), max(acc_count_beta / iter)
        ))
      }
    }

    # finalize posterior mean of pi
    ##n_saved <- max(1, iter_save)
    out_pi_mean <- out_pi_mean / max(1, nsave) # average over saved iterations

    ret <- list(
      Beta_samples = out_beta,
      nu = out_nu,
      Sigma = out_Sigma,
      z_samples = out_z,
      pi_ik = out_pi_ik,
      pi_mean = out_pi_mean,
      estimate_nu = estimate_nu,
      burnin = burnin, niter = niter, thin = thin,
      loglik = logliks, loglik_individual = logliks_individual
    )
    class(ret) <- c("moewishart.bayes")
  }

  if (method == "em") {
    maxit <- niter
    ret <- moewishart.em(
      S_list, X, K, maxit, tol, verbose,
      init, estimate_nu, init_nu, ridge
    )
    ret$estimate_nu <- estimate_nu
    class(ret) <- c("moewishart.em")
  }

  return(ret)
}

# -- Internal help functions for full Bayesian--

# helper: compute pi_ik matrix (n x K) given Beta
softmax_rows <- function(L) {
  # L: n x K matrix of linear predictors
  m <- apply(L, 1, max)
  Ls <- L - m
  S <- exp(Ls)
  row_sums <- rowSums(S)
  S / row_sums
}
compute_pi_ik <- function(X, Beta) {
  L <- X %*% Beta # n x K
  # numeric stability: subtract row max
  rm <- apply(L, 1, max)
  Ls <- L - rm
  expL <- exp(Ls)
  expL / rowSums(expL)
}

# -- Internal function with EM algorithm for the Wishart mixture-of-experts model--

moewishart.em <- function(S_list, X, K, maxit = 200, tol = 1e-6, verbose = TRUE,
                          init = NULL, estimate_nu = FALSE, init_nu = NULL, ridge = 1e-8) {
  n <- length(S_list)
  p <- nrow(S_list[[1]])
  q <- ncol(X)

  # --- CORRECTION 1: BREAK SYMMETRY IN INITIALIZATION ---

  if (is.null(init) || is.null(init$beta)) {
    beta_vec <- rep(0, q * (K - 1))
  } else {
    beta_vec <- as.numeric(init$beta)
  }

  if (is.null(init) || is.null(init$Sigma)) {
    pooled <- Reduce("+", S_list) / n
    if (is.null(init_nu)) init_nu <- rep(p + 5, K)
    # Add slight jitter to initial Sigmas to ensure they aren't identical
    Sigma_list <- lapply(1:K, function(k) {
      (pooled / init_nu[k]) * (1 + runif(1, -0.1, 0.1))
    })
  } else {
    Sigma_list <- init$Sigma
  }

  if (is.null(init_nu)) init_nu <- rep(p + 5, K)
  nu_vec <- init_nu

  loglik_hist <- numeric(maxit)
  small_eps <- ridge

  # precompute log|S_i| for all i
  logdetS <- sapply(S_list, function(S) as.numeric(determinant(S, logarithm = TRUE)$modulus))

  # EM loop
  for (iter in 1:maxit) {
    # --- E-step ---
    pi_mat <- gating_probs_from_beta(beta_vec, q, K, X)
    log_f_mat <- matrix(NA, n, K)

    for (k in 1:K) {
      for (i in 1:n) {
        # Use precomputed logdetS to speed up
        log_f_mat[i, k] <- dWishart(S_list[[i]], nu_vec[k], Sigma_list[[k]],
          logdetS[i],
          logarithm = TRUE
        )
      }
    }

    log_post <- log(pi_mat + 1e-300) + log_f_mat
    row_max <- apply(log_post, 1, max)
    exp_shifted <- exp(log_post - row_max)
    row_sums <- rowSums(exp_shifted)
    gamma <- exp_shifted / row_sums

    loglik <- sum(row_max + log(row_sums))
    loglik_hist[iter] <- loglik

    if (verbose) cat(sprintf("Iter %3d  loglik = %.6f\n", iter, loglik))

    if (iter > 1 && abs(loglik - loglik_hist[iter - 1]) < tol) {
      if (verbose) cat("Converged by loglik tolerance.\n")
      break
    }

    # --- M-step ---
    an <- compute_A_nk(gamma, K, p, n, S_list)
    A_list <- an$A
    n_k <- an$n_k

    # Avoid division by zero if a cluster dies
    valid_k <- n_k > 1e-6

    for (k in 1:K) {
      if (valid_k[k]) {
        if (estimate_nu) {
          mean_logS_k <- sum(gamma[, k] * logdetS) / n_k[k]
          nu_old <- nu_vec[k]
          nu_vec[k] <- solve_nu(nu_old, A_list[[k]], n_k[k], mean_logS_k, p)
        }

        # Update Sigma using the NEW nu
        Sigma_new <- A_list[[k]] / (nu_vec[k] * n_k[k])
        Sigma_new <- Sigma_new + diag(small_eps, p)
        Sigma_list[[k]] <- (Sigma_new + t(Sigma_new)) / 2
      } else {
        # Reset empty component (optional strategy)
        if (verbose) cat(sprintf("Resetting component %d\n", k))
        nu_vec[k] <- p + 5
        Sigma_list[[k]] <- diag(1, p)
      }
    }

    # Update beta
    opt <- optim(beta_vec,
      fn = neg_wt_multinom, gr = grad_wt_multinom,
      q, K, X, # Further arguments to be passed to fn and gr
      method = "BFGS", gamma = gamma, control = list(maxit = 50)
    )
    beta_vec <- opt$par
  }

  Beta_mat <- matrix(0, q, K)
  if (K > 1) Beta_mat[, 1:(K - 1)] <- matrix(beta_vec, q, K - 1)
  colnames(Beta_mat) <- paste0("comp", 1:K)

  list(
    K = K, p = p, q = q, n = n,
    Beta = Beta_mat, Sigma = Sigma_list, nu = nu_vec,
    gamma = gamma, loglik = loglik_hist[1:iter], iter = iter
  )
}

# -- Internal help functions for EM algorithm--

# helper: log multivariate gamma derivative (psi_p) and its derivative
psi_p <- function(a, p) sum(digamma(a + (1 - (1:p)) / 2))

trigamma_p <- function(a, p) sum(trigamma(a + (1 - (1:p)) / 2))

compute_A_nk <- function(gamma, K, p, n, S_list) {
  n_k <- colSums(gamma)
  A_list <- vector("list", K)
  for (k in 1:K) {
    A <- matrix(0, p, p)
    # Vectorizing this sum is hard with list of matrices, loop is okay
    for (i in 1:n) A <- A + gamma[i, k] * S_list[[i]]
    A_list[[k]] <- A
  }
  list(A = A_list, n_k = n_k)
}
gating_probs_from_beta <- function(beta_vec, q, K, X) {
  B <- matrix(0, q, K)
  if (K > 1) {
    B[, 1:(K - 1)] <- matrix(beta_vec, q, K - 1)
  }
  eta <- X %*% B
  eta_max <- apply(eta, 1, max)
  exp_eta <- exp(eta - eta_max)
  denom <- rowSums(exp_eta)
  exp_eta / denom
}

neg_wt_multinom <- function(beta_vec, gamma, q, K, X) {
  pi_mat <- gating_probs_from_beta(beta_vec, q, K, X)
  # prevent log(0)
  ll <- sum(gamma * log(pi_mat + 1e-300))
  -ll
}

grad_wt_multinom <- function(beta_vec, gamma, q, K, X) {
  B <- matrix(0, q, K)
  if (K > 1) B[, 1:(K - 1)] <- matrix(beta_vec, q, K - 1)
  pi_mat <- gating_probs_from_beta(beta_vec, q, K, X)
  grad_mat <- matrix(0, q, K - 1)
  for (k in 1:(K - 1)) {
    # gradient of NEGATIVE likelihood is sum((pi - gamma) * X)
    diff <- pi_mat[, k] - gamma[, k]
    grad_mat[, k] <- colSums(X * diff)
  }
  as.numeric(grad_mat)
}

# --- CORRECTION 2: FIX NEWTON-RAPHSON DERIVATIVE ---
solve_nu <- function(init_nu, A_k, n_k, mean_logS_k, p, maxit = 20) {
  nu <- max(init_nu, p + 1.001)

  # OPTIMIZATION: Compute logdet(A) once outside the loop
  # A_k is constant during the nu update
  logdetA <- as.numeric(determinant(A_k, logarithm = TRUE)$modulus)

  for (it in 1:maxit) {
    # Use algebraic property: log|A / (n*nu)| = log|A| - p*log(n*nu)
    logdetSigma <- logdetA - p * log(n_k * nu)

    # Objective function
    lhs <- psi_p(nu / 2, p)
    rhs <- mean_logS_k - logdetSigma - p * log(2)
    fval <- lhs - rhs

    if (abs(fval) < 1e-6) break

    # Derivative
    df <- 0.5 * trigamma_p(nu / 2, p) - (p / nu)

    # Newton update
    nu_new <- nu - fval / df

    # Boundary check
    if (!is.finite(nu_new) || nu_new <= p + 1e-6) {
      nu_new <- (nu + (p + 1.0)) / 2
    }
    if (abs(nu_new - nu) < 1e-6) {
      nu <- nu_new
      break
    }
    nu <- nu_new
  }
  max(nu, p + 1.001)
}
