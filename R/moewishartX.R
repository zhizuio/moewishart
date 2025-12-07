#' @title Algorithm for Wishart mixture(-of-experts) model with covariates
#'
#' @description
#' TBA
#'
#' @name moewishartX
#'
#' @importFrom utils combn
#' @importFrom stats optim
#'
#' @param S_list TBA
#' @param X TBA
#' @param K TBA
#' @param niter TBA
#' @param burnin TBA
#' @param method one of \code{c("bayes", "em")}
#' @param thin TBA
#' @param nu0 TBA
#' @param Psi0 TBA
#' @param init_nu TBA
#' @param estimate_nu TBA
#' @param nu_prior_a TBA
#' @param nu_prior_b TBA
#' @param mh_sigma TBA
#' @param mh_beta TBA
#' @param sigma_beta TBA
#' @param init TBA
#' @param tol TBA
#' @param ridge TBA
#' @param verbose TBA
#'
#'
#' @return An object of a list including ...
#'
#'
#' @examples
#'
#' # simulate data
#' set.seed(123)
#' n <- 200 # subjects
#' p <- 10
#'
#' @export
moewishartX <- function(S_list,
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
    nsave <- floor((niter - burnin) / thin)
    if (nsave < 1) nsave <- 1
    out_beta <- array(NA, dim = c(nsave, q, K)) # store Beta at saves
    out_nu <- matrix(NA, nrow = nsave, ncol = K)
    out_Sigma <- vector("list", nsave)
    out_z <- matrix(NA, nrow = nsave, ncol = n)
    out_pi_mean <- array(0, dim = c(n, K)) # accumulate posterior mean of pi_ik
    logliks <- numeric(niter)
    iter_save <- 0

    # pre-alloc
    logpost <- matrix(0, n, K)
    acc_count_nu <- numeric(K) # record acceptance counts of MH for nu
    acc_count_beta <- numeric(K - 1) # record acceptance counts of MH for beta

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
      free_idx_cols <- 1:(K - 1)
      Beta_free <- as.vector(Beta[, free_idx_cols]) # length q*(K-1)
      prop <- Beta_free + rnorm(length(Beta_free), 0, mh_beta)
      Beta_prop <- Beta
      Beta_prop[, free_idx_cols] <- matrix(prop, nrow = q, ncol = K - 1)
      Beta_prop[, K] <- 0
      # compute new pi_ik (n x K)
      pi_prop <- compute_pi_ik(X, Beta_prop)
      # log prior for Beta (Gaussian iid)
      lp_prior_old <- -0.5 * sum(Beta_free^2) / (sigma_beta^2)
      lp_prior_new <- -0.5 * sum(Beta_prop^2) / (sigma_beta^2)
      # log-likelihood of labels given gating: sum_i log pi_i,z[i] (only depends on z)
      ll_old <- sum(log(pi_ik[cbind(1:n, z)] + 1e-300))
      ll_new <- sum(log(pi_prop[cbind(1:n, z)] + 1e-300))
      log_accept <- (ll_new + lp_prior_new) - (ll_old + lp_prior_old)
      if (log(runif(1)) < log_accept) {
        Beta <- Beta_prop
        pi_ik <- pi_prop
        # optionally track acceptance
        acc_count_beta <- acc_count_beta + 1
      }

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
          prop_log <- rnorm(1, log(curr_nu), mh_sigma)
          prop_nu <- exp(prop_log)
          if (prop_nu > p - 1 + 1e-8) {
            idx <- which(z == k)
            if (length(idx) > 0) {
              sum_log_det_S_k <- sum(log_det_S[idx])
              Sig <- Sigma_k[, , k]
              chol_Sig <- tryCatch(chol(Sig), error = function(e) chol(Sig + diag(1e-8, p)))
              log_det_Sig <- 2 * sum(log(diag(chol_Sig)))
              Sig_inv <- chol2inv(chol_Sig)
              S_sum_vec <- colSums(S_mat[idx, , drop = FALSE])
              sum_tr_val <- sum(as.vector(Sig_inv) * S_sum_vec)
              calc_ll_nu <- function(val_nu) {
                term1 <- (val_nu - p - 1) / 2 * sum_log_det_S_k
                term2 <- -0.5 * sum_tr_val
                term3 <- -length(idx) * ((val_nu * p / 2) * log(2) + (val_nu / 2) * log_det_Sig)
                term4 <- -length(idx) * lmvgamma(val_nu / 2, p)
                term1 + term2 + term3 + term4
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

      # --- Save samples after burnin and thinning ---
      if (iter > burnin && ((iter - burnin) %% thin == 0)) {
        iter_save <- iter_save + 1
        if (iter_save <= nsave) {
          out_beta[iter_save, , ] <- Beta
          out_nu[iter_save, ] <- nu_k
          out_Sigma[[iter_save]] <- Sigma_k
          out_z[iter_save, ] <- z
        }
        out_pi_mean <- out_pi_mean + pi_ik
      }

      if (verbose && (iter %% 500 == 0 || iter == 1)) {
        # elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        # rate <- iter / elapsed
        cat(sprintf(
          "Iter %4d | LL=%.1f | acc_rate_nu=[%.3f %.3f] | acc_rate_beta=%.3f\n",
          iter, logliks[iter], min(acc_count_nu / iter), max(acc_count_nu / iter),
          acc_count_beta / iter
        ))
      }
    }

    # finalize posterior mean of pi
    n_saved <- max(1, iter_save)
    out_pi_mean <- out_pi_mean / max(1, nsave) # average over saved iterations

    ret <- list(
      Beta_samples = out_beta,
      nu_samples = out_nu,
      Sigma_samples = out_Sigma,
      z_samples = out_z,
      pi_mean = out_pi_mean,
      loglik = logliks
    )
  }

  if (method == "em") {
    maxit <- niter
    ret <- moewishartX.em(
      S_list, X, K, maxit, tol, verbose,
      init, estimate_nu, init_nu, ridge
    )
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

moewishartX.em <- function(S_list, X, K, maxit = 200, tol = 1e-6, verbose = TRUE,
                           init = NULL, estimate_nu = FALSE, init_nu = NULL, ridge = 1e-8) {
  n <- length(S_list)
  p <- nrow(S_list[[1]])
  q <- ncol(X)

  # --- CORRECTION 1: BREAK SYMMETRY IN INITIALIZATION ---
  if (is.null(init)) {
    # Initialize with random noise so clusters can diverge
    gamma <- matrix(runif(n * K), n, K)
    gamma <- gamma / rowSums(gamma)

    if (is.null(init_nu)) init_nu <- rep(p + 5, K)
  } else {
    gamma <- init$gamma
    if (is.null(init_nu) && !is.null(init$nu)) init_nu <- init$nu
  }

  pars <- list()
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
