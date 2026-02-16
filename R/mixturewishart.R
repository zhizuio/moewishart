#' @title EM/Bayesian estimation for Wishart mixture model
#'
#' @description
#' Fit finite mixtures of Wishart-distributed SPD matrices using either a
#' Bayesian sampler or the EM algorithm. The input \code{S_list} is a list
#' of \eqn{p \times p} SPD matrices. Under component \eqn{k},
#' \eqn{S_i \mid z_i=k \sim W_p(\nu_k, \Sigma_k)} with degrees of freedom
#' \eqn{\nu_k} and SPD scale matrix \eqn{\Sigma_k}. Mixture weights
#' \eqn{\pi_k} sum to \eqn{1}.
#'
#' @name mixturewishart
#'
#' @importFrom utils combn
#' @importFrom stats kmeans
#'
#' @param S_list List of length \eqn{n} of SPD matrices, each \eqn{p \times p}.
#'   These are the observed matrices modeled by a mixture of Wisharts.
#' @param K Integer. Number of mixture components.
#' @param niter Integer. Total iterations. Bayesian mode: total MCMC
#'   iterations (including burn-in). EM mode: maximum EM iterations
#'   (alias to \code{maxiter}).
#' @param burnin Integer. Number of burn-in iterations (Bayesian mode).
#' @param method Character; one of \code{c("bayes","em")}. Selects sampler
#'   or optimizer.
#' @param thin Integer. Thinning interval for saving draws (Bayesian).
#' @param alpha Numeric vector length \eqn{K} (Dirichlet prior on
#'   \eqn{\pi}) or \code{NULL} to default to \code{rep(1, K)} (Bayesian).
#' @param nu0 Numeric. Inverse-Wishart prior df for \eqn{\Sigma_k}
#'   (Bayesian). Default: \eqn{p + 2}.
#' @param Psi0 Numeric \eqn{p \times p} SPD matrix. Inverse-Wishart prior
#'   scale for \eqn{\Sigma_k} (Bayesian). Default: \code{diag(p)}.
#' @param init_pi Optional numeric vector length \eqn{K} summing to
#'   \eqn{1}. EM initialization for mixture weights. If \code{NULL},
#'   random or data-driven initialization is used.
#' @param init_nu Optional numeric vector length \eqn{K} of initial
#'   degrees of freedom. Used in both modes.
#' @param init_Sigma Optional list of \eqn{K} SPD matrices (each
#'   \eqn{p \times p}). EM initialization for \eqn{\Sigma_k}.
#' @param marginal.z Logical. If \code{TRUE}, integrates out \eqn{\pi}
#'   when sampling \eqn{z} (collapsed step) in Bayesian mode. If
#'   \code{FALSE}, samples \eqn{z} conditional on current \eqn{\pi}.
#' @param estimate_nu Logical. If \code{TRUE}, estimate/update
#'   \eqn{\nu_k} (MH in Bayesian mode; Newton/EM in EM). If
#'   \code{FALSE}, \eqn{\nu_k} are fixed.
#' @param nu_prior_a Numeric. Prior hyperparameter \eqn{a} for
#'   \eqn{\nu_k} (Bayesian), used when \code{estimate_nu = TRUE}.
#' @param nu_prior_b Numeric. Prior hyperparameter \eqn{b} for
#'   \eqn{\nu_k} (Bayesian), used when \code{estimate_nu = TRUE}.
#' @param mh_sigma Numeric scalar or length-\eqn{K} vector. Proposal sd
#'   for MH updates on \eqn{\log(\nu_k)} (Bayesian, when estimating
#'   \eqn{\nu}).
#' @param n_restarts Integer. Number of random restarts for EM. Ignored
#'   in Bayesian mode.
#' @param restart_iters Integer. Number of short EM iterations per
#'   restart used to select a good initialization. Ignored in Bayesian
#'   mode.
#' @param tol Numeric. Convergence tolerance on absolute change of
#'   log-likelihood (EM), also used internally elsewhere.
#' @param verbose Logical. If \code{TRUE}, print progress information.
#'
#'
#' @details
#' Mixture mixture model:
#'  \eqn{p(S_i) = \sum_{k=1}^K \pi_k \, f_W(S_i \mid \nu_k, \Sigma_k)}.
#'
#' Algorithms:
#' \enumerate{
#'   \item \code{method = "bayes"}: Samples latent labels \eqn{z}, weights
#'         \eqn{\pi}, component scales \eqn{\Sigma_k}, and optionally
#'         \eqn{\nu_k}. Uses a Dirichlet prior for \eqn{\pi}, inverse-
#'         Wishart prior for \eqn{\Sigma_k}, and a prior on \eqn{\nu_k}
#'         when \code{estimate_nu = TRUE}. Degrees-of-freedom are updated
#'         via MH on \eqn{\log(\nu_k)} with proposal sd \code{mh_sigma}.
#'         Can integrate out \eqn{\pi} when sampling \eqn{z} if
#'         \code{marginal.z = TRUE}.
#'   \item \code{method = "em"}: Maximizes the observed-data log-
#'         likelihood via EM. The E-step computes responsibilities via
#'         Wishart log-densities. The M-step updates \eqn{\pi_k} and
#'         \eqn{\Sigma_k}; optionally updates \eqn{\nu_k} when
#'         \code{estimate_nu = TRUE}. Supports multiple random restarts.
#' }
#'
#' Note that
#' (i) All matrices in \code{S_list} must be SPD. Small ridge terms may be
#' added internally for stability, and
#' (ii) Multiple EM restarts are recommended for robustness on difficult datasets.
#'
#'
#' @return A list whose structure depends on \code{method}:
#' \itemize{
#'   \item For \code{method = "bayes"}:
#'     \itemize{
#'       \item \code{pi_ik}: array (\code{nsave} x \code{n} x \code{K}),
#'             saved per-observation weights.
#'       \item \code{pi}: matrix (\code{nsave} x \code{K}), saved mixture
#'             proportions.
#'       \item \code{nu}: matrix (\code{nsave} x \code{K}), saved degrees-
#'             of-freedom.
#'       \item \code{Sigma}: list of length \code{nsave}; each is an
#'             array (\eqn{p \times p \times K}) of \eqn{\Sigma_k} draws.
#'       \item \code{z}: matrix (\code{nsave} x \code{n}), saved
#'             allocations.
#'       \item \code{sigma_posterior_mean}: array (\eqn{p \times p \times
#'             K}), posterior mean of \eqn{\Sigma_k}.
#'       \item \code{loglik}: numeric vector (length \code{niter}), log-
#'             likelihood trace.
#'       \item \code{loglik_individual}: matrix (\code{niter} x
#'             \code{n}), per-observation log-likelihood.
#'     }
#'   \item For \code{method = "em"}:
#'     \itemize{
#'       \item \code{pi}: numeric vector length \eqn{K}, mixture
#'             proportions.
#'       \item \code{Sigma}: list length \code{K}, each a \eqn{p \times p}
#'             SPD matrix.
#'       \item \code{nu}: numeric vector length \code{K}, degrees-of-
#'             freedom.
#'       \item \code{tau}: matrix (\eqn{n \times K}), responsibilities.
#'       \item \code{loglik}: numeric vector, log-likelihood per EM
#'             iteration.
#'       \item \code{iterations}: integer, number of EM iterations
#'             performed.
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
#' dat <- simData(n, p,
#'   K = 3,
#'   pis = c(0.35, 0.40, 0.25),
#'   nus = c(8, 16, 3)
#' )
#'
#' set.seed(123)
#' fit <- mixturewishart(
#'   dat$S,
#'   K = 3,
#'   mh_sigma = c(0.2, 0.1, 0.1), # tune this for MH acceptance 20-40%
#'   niter = 100, burnin = 50
#' )
#'
#' # Posterior means for degrees of freedom of Wishart distributions:
#' nu_mcmc <- fit$nu[-c(1:fit$burnin), ]
#' colMeans(nu_mcmc)
#'
#' @export
mixturewishart <- function(S_list,
                           K,
                           niter = 3000,
                           burnin = 1000,
                           method = "bayes",
                           thin = 1,
                           alpha = NULL,
                           nu0 = NULL,
                           Psi0 = NULL,
                           init_pi = NULL,
                           init_nu = NULL,
                           init_Sigma = NULL,
                           marginal.z = TRUE,
                           estimate_nu = TRUE,
                           nu_prior_a = 2,
                           nu_prior_b = 0.1,
                           mh_sigma = 1,
                           n_restarts = 3,
                           restart_iters = 20,
                           tol = 1e-6,
                           verbose = TRUE) {
  if (!method %in% c("bayes", "em")) {
    stop("Argument 'method' must be either 'bayes' or 'em'!")
  }

  # TODO: remove redundant code for the common calculations between full Bayesian and EM algorithm
  if (method == "bayes") {
    # -- 1. Pre-processing and Pre-allocation --
    n <- length(S_list)
    p <- nrow(S_list[[1]])

    if (length(mh_sigma) != K && length(mh_sigma) != 1) {
      stop("Argument 'mh_sigma' must have length 1 or K!")
    } else if (length(mh_sigma) == 1) {
      mh_sigma <- rep(mh_sigma, K)
    }

    if (!is.null(alpha)) {
      if (length(alpha) != K) {
        warning("Length of alpha (", length(alpha), ") != K (", K, "). Recycling/triming alpha to length K.")
        alpha <- rep(alpha, length.out = K)
      }
    } else {
      alpha <- rep(1, K)
    }

    # Defaults
    if (is.null(nu0)) nu0 <- p + 2
    if (is.null(Psi0)) Psi0 <- diag(p)
    if (is.null(init_nu)) init_nu <- rep(p + 2, K)

    # OPTIMIZATION: Vectorize Data
    # Flatten each p x p matrix into a row of length p^2
    # This allows fast summation and fast trace calculation
    S_mat <- t(sapply(S_list, as.vector)) # Dimension: n x (p*p)

    # OPTIMIZATION: Precompute log determinants of data
    # This part of the density never changes
    log_det_S <- sapply(S_list, function(x) determinant(x, logarithm = TRUE)$modulus)

    # Initialize Parameters
    # Use vectorized data for kmeans
    km <- kmeans(S_mat, centers = K, nstart = 5)
    z <- km$cluster

    if (is.null(init_pi)) {
      pi_k <- table(factor(z, levels = 1:K)) / n
    } else {
      if (length(init_pi) != K || sum(init_pi) != 1) {
        stop("Please specify correct 'init_pi'!")
      }
      pi_k <- init_pi
    }
    Sigma_k <- array(0, c(p, p, K))
    nu_k <- init_nu

    # Initial Sigma
    for (k in 1:K) {
      idx <- which(z == k)
      if (length(idx) == 0) {
        Sigma_k[, , k] <- Psi0 / (nu0 - p - 1)
      } else {
        # OPTIMIZATION: Fast matrix sum using colSums on vectorized data
        S_sum_vec <- colSums(S_mat[idx, , drop = FALSE])
        S_sum <- matrix(S_sum_vec, p, p)
        Sigma_k[, , k] <- (Psi0 + S_sum) / (nu0 + length(idx) * nu_k[k] - p - 1)
      }
    }

    # Storage
    ## nsave <- floor((niter - burnin) / thin)
    nsave <- floor(niter / thin)
    if (nsave < 1) nsave <- 1
    out_pi_ik <- array(NA, dim = c(nsave, n, K))
    out_pi <- matrix(NA, nrow = nsave, ncol = K)
    out_nu <- matrix(NA, nrow = nsave, ncol = K)
    out_Sigma <- vector("list", nsave)
    out_z <- matrix(NA, nrow = nsave, ncol = n)
    logliks <- numeric(niter)
    logliks_individual <- matrix(NA, nrow = niter, ncol = n)
    iter_save <- 0

    # Pre-allocate reusable vectors
    logpost <- matrix(0, n, K)
    n_k <- as.numeric(table(factor(z, levels = 1:K)))
    acc_count <- numeric(K) # record acceptance counts of MH for nu

    # Start Timer
    # start_time <- Sys.time()

    for (iter in 1:niter) {
      # --- Step 1: Update Labels z (The Heavy Lifting) ---

      for (k in 1:K) {
        # OPTIMIZATION: Invert Sigma ONLY ONCE per cluster
        Sig <- Sigma_k[, , k]

        # Cholesky is faster and more stable for determinant/inverse
        chol_Sig <- tryCatch(chol(Sig), error = function(e) chol(Sig + diag(1e-6, p)))
        log_det_Sig <- 2 * sum(log(diag(chol_Sig))) # log|Sigma|
        Sig_inv <- chol2inv(chol_Sig) # Sigma^-1

        # OPTIMIZATION: Vectorized Trace calculation
        # Tr(Sigma^-1 * S_i) is the dot product of vec(Sigma^-1) and vec(S_i)
        # We calculate this for ALL i at once via matrix multiplication
        # S_mat is (n x p^2), as.vector(Sig_inv) is (p^2 x 1) -> Result (n x 1)
        tr_val <- S_mat %*% as.vector(Sig_inv)

        # Calculate log-density for all n points
        # term1: (nu - p - 1)/2 * log|S|
        # term2: -0.5 * tr(Sig^-1 S)
        # term3: Normalizing constants involving nu and log|Sig|

        nu <- nu_k[k]

        term1 <- (nu - p - 1) / 2 * log_det_S
        term2 <- -0.5 * tr_val
        term3 <- -(nu * p / 2) * log(2) - (nu / 2) * log_det_Sig
        term4 <- -lmvgamma(nu / 2, p)

        # logpost[, k] <- log(pi_k[k] + 1e-300) + term1 + term2 + term3 + term4
        logpost[, k] <- term1 + term2 + term3 + term4
        if (marginal.z) {
          logpost[, k] <- logpost[, k] + log(alpha[k] + n_k[k] - as.numeric(z == k) + 1e-300)
        } else {
          logpost[, k] <- logpost[, k] + log(pi_k[k] + 1e-300)
        }
      }

      # Sample z
      # Vectorized sampling is hard in base R, looping sample.int is okay
      # but we can optimize the probability normalization
      # Using pure R loop for sampling is usually fast enough compared to the math above
      pi_ik <- matrix(NA, nrow = n, ncol = K)
      for (i in 1:n) {
        lp <- logpost[i, ]
        lp <- lp - max(lp)
        prob <- exp(lp)
        z[i] <- sample.int(K, 1, prob = prob)
        pi_ik[i, ] <- prob / sum(prob)
      }

      # --- Step 2: Update Weights pi ---
      n_k <- as.numeric(table(factor(z, levels = 1:K)))
      pi_k <- as.numeric(rdirichlet(1, alpha + n_k))
      # pi_k <- c(0.35, 0.40, 0.25)

      # --- Step 3: Update Sigma_k ---
      for (k in 1:K) {
        idx <- which(z == k)
        nk <- length(idx)

        if (nk == 0) {
          Sigma_k[, , k] <- sampleIW(nu0, solve(Psi0))
        } else {
          # OPTIMIZATION: Fast Sum
          S_sum_vec <- colSums(S_mat[idx, , drop = FALSE])
          S_sum <- matrix(S_sum_vec, p, p)

          nu_post <- nu0 + nk * nu_k[k]
          Psi_post <- solve(Psi0) + S_sum

          # Invert Psi_post once for sampling
          # Use tryCatch for numerical stability
          Psi_post_inv <- tryCatch(solve(Psi_post), error = function(e) solve(Psi_post + diag(1e-6, p)))

          # browser() ##TODO: check bugs in updating Sigma_k
          Sigma_k[, , k] <- sampleIW(nu_post, Psi_post_inv)
          # if (k == 1 )
          #   Sigma_k[, , k] <- matrix(c(0.5,0.2,0.2,0.7), nrow = 2)
          # if (k == 2 )
          #   Sigma_k[, , k] <- matrix(c(2,0.6,0.6,1.5), nrow = 2)
          # if (k == 3 )
          #   Sigma_k[, , k] <- matrix(c(4,0.2,0.2,3), nrow = 2)
        }
      }

      # --- Step 4: Update nu (MH) ---
      if (estimate_nu) {
        for (k in 1:K) {
          curr_nu <- nu_k[k]
          prop_log <- rnorm(1, log(curr_nu), mh_sigma[k])
          # prop_log <- log(curr_nu) + rnorm(1, 0, mh_sigma) # random-walk MH
          prop_nu <- exp(prop_log)

          if (prop_nu > p - 1 + 1e-6) {
            # We need the Likelihood sum for this cluster
            # Reuse the data stats we already know
            idx <- which(z == k)
            if (length(idx) > 0) {
              # Re-calculate only necessary parts
              # We need log|S| sum and tr(Sig^-1 S) sum

              # Grab the specific rows and sum them up for efficiency
              sum_log_det_S_k <- sum(log_det_S[idx])

              # We already have Sig_inv from Step 1?
              # No, Step 3 updated Sigma. We must re-invert current Sigma.
              Sig <- Sigma_k[, , k]
              chol_Sig <- tryCatch(chol(Sig), error = function(e) chol(Sig + diag(1e-6, p)))
              log_det_Sig <- 2 * sum(log(diag(chol_Sig)))
              # Sig_inv <- chol2inv(chol_Sig)

              # # Sum of traces = Tr(Sig^-1 * Sum(S))
              # # We can calculate Sum(S) fast
              # S_sum_vec <- colSums(S_mat[idx, , drop = FALSE])
              # sum_tr_val <- sum(as.vector(Sig_inv) * S_sum_vec)

              # Define a mini function for log-lik given Sufficient Stats
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

            # Priors + Jacobian
            lp_old <- (nu_prior_a - 1) * log(curr_nu) - nu_prior_b * curr_nu + log(curr_nu)
            lp_new <- (nu_prior_a - 1) * log(prop_nu) - nu_prior_b * prop_nu + log(prop_nu)

            if (log(runif(1)) < (ll_new + lp_new) - (ll_old + lp_old)) {
              nu_k[k] <- prop_nu
              acc_count[k] <- acc_count[k] + 1
            }
          }
        }
      }
      # browser() ##TODO: check bugs in updating nu_k
      # nu_k <- c(8, 12, 3)

      # --- Calculate LogLik for history (fast approximation using Step 1 data) ---
      # We actually calculated logpost at the start of the loop (using old params).
      # We can use that for the record, or re-calc.
      # Using the one from Step 1 is "Lag-1" loglik but much faster.
      # For strictness, let's recalculate using log-sum-exp on logpost:
      # Note: logpost updated in Step 1 corresponds to Z sampling.

      # Fast LogSumExp on rows
      max_l <- apply(logpost, 1, max)
      row_sums <- exp(logpost - max_l)
      logliks[iter] <- sum(max_l + log(rowSums(row_sums)))
      logliks_individual[iter, ] <- max_l + log(rowSums(row_sums))

      # --- Save ---
      ## if (iter > burnin && ((iter - burnin) %% thin == 0)) {
      if (iter %% thin == 0) {
        iter_save <- iter_save + 1
        if (iter_save <= nsave) {
          out_pi_ik[iter_save, , ] <- pi_ik
          out_pi[iter_save, ] <- pi_k
          out_nu[iter_save, ] <- nu_k
          out_Sigma[[iter_save]] <- Sigma_k
          out_z[iter_save, ] <- z
        }
      }

      if (verbose && (iter %% 500 == 0 || iter == 1)) {
        # Calculate speed
        # elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        # rate <- iter / elapsed
        cat(sprintf(
          "Iter %4d | LL=%.1f | acc_rate_nu=%.3f\n",
          iter, logliks[iter], acc_count / iter # min(acc_count / iter), max(acc_count / iter)
        ))
      }
    }

    sigma_posterior_mean <- Reduce("+", out_Sigma) / length(out_Sigma)

    ret <- list(
      pi_ik = out_pi_ik,
      pi = out_pi, nu = out_nu, Sigma = out_Sigma, z = out_z,
      sigma_posterior_mean = sigma_posterior_mean,
      estimate_nu = estimate_nu,
      burnin = burnin, niter = niter, thin = thin,
      loglik = logliks,
      loglik_individual = logliks_individual
    )

    class(ret) <- c("mixturewishart.bayes")
  }

  if (method == "em") {
    maxiter <- niter
    ret <- mixturewishart.em(
      S_list, K,
      n_restarts, restart_iters,
      init_pi, init_Sigma, init_nu,
      maxiter, tol, estimate_nu, verbose
    )

    ret$estimate_nu <- estimate_nu
    class(ret) <- c("mixturewishart.em")
  }

  return(ret)
}

# -- Internal function with EM algorithm for the Wishart mixture model--

mixturewishart.em <- function(S_list, K,
                              n_restarts = 3,
                              restart_iters = 20,
                              init_pi = NULL,
                              init_Sigma = NULL,
                              init_nu = NULL,
                              maxiter = 200,
                              tol = 1e-6,
                              estimate_nu = TRUE,
                              verbose = TRUE) {
  n <- length(S_list)
  p <- ncol(S_list[[1]])

  # Pre-compute log|S_i| as it is constant throughout EM
  logdetS <- sapply(S_list, function(S) {
    as.numeric(determinant(S + diag(p) * 1e-12, logarithm = TRUE)$modulus)
  })

  # --- 1. Initialization (CORRECTED) ---

  # =========================================================================
  # PHASE 1: Initialization Selection (Multiple Restarts)
  # =========================================================================

  # We only run restarts if the user did NOT provide specific starting Sigma/pi
  if (is.null(init_Sigma) && is.null(init_pi) && is.null(init_nu)) {
    if (verbose) cat("Running", n_restarts, "initialization restarts...\n")

    best_start_loglik <- -Inf
    best_params <- list()

    r <- 1
    repeat {
      ## for (r in 1:n_restarts) {

      # A. Random Initialization for this attempt
      # Default nu
      curr_nu <- if (is.null(init_nu)) rep_len(c(p + 1, p + 5), length.out = K) else init_nu
      tamtam <- rgamma(K, shape = 1)
      curr_pi <- tamtam / sum(tamtam)

      # Random centers
      rand_indices <- sample(seq_len(n), K)
      curr_Sigma <- lapply(seq_len(K), function(k) {
        idx <- rand_indices[k]
        S_list[[idx]] / curr_nu[k]
      })

      # B. Run "Small EM" (Short Loop)
      curr_loglik <- -Inf

      for (small_it in 1:restart_iters) {
        # E-Step (Simplified)
        logdens <- matrix(NA, n, K)
        for (k in seq_len(K)) {
          for (i in seq_len(n)) {
            # if(is.null(curr_Sigma[[k]])) browser()
            logdens[i, k] <- dWishart(S_list[[i]], curr_nu[k], curr_Sigma[[k]],
              logarithm = TRUE
            )
          }
          logdens[, k] <- logdens[, k] + log(curr_pi[k])
        }
        maxlog <- apply(logdens, 1, max)
        denom <- maxlog + log(rowSums(exp(logdens - maxlog)))
        logtau <- sweep(logdens, 1, denom, FUN = "-")
        tau <- exp(logtau)
        curr_loglik <- sum(denom)

        # M-Step (Simplified - update Pi and Sigma only, keep Nu fixed for stability in init)
        N_k <- colSums(tau)
        curr_pi <- N_k / n

        # Check for collapse
        if (any(N_k < 1e-6)) {
          curr_loglik <- -Inf
          break
        }

        for (k in seq_len(K)) {
          Ssum <- matrix(0, p, p)
          for (i in seq_len(n)) if (tau[i, k] > 1e-10) Ssum <- Ssum + tau[i, k] * S_list[[i]]
          curr_Sigma[[k]] <- Ssum / (N_k[k] * curr_nu[k])
        }
      }

      # C. Compare and Store
      if (verbose) cat(sprintf("  -> Restart %d: Loglik = %.2f\n", r, curr_loglik))

      if (curr_loglik > best_start_loglik) {
        best_start_loglik <- curr_loglik
        best_params <- list(pi = curr_pi, Sigma = curr_Sigma, nu = curr_nu)
      }

      if (r >= n_restarts && curr_loglik > -Inf) break
      r <- r + 1
    }

    # Set the main loop variables to the winner
    pi_k <- best_params$pi
    Sigma_k <- best_params$Sigma
    nu_k <- best_params$nu
  } else {
    if (is.null(init_Sigma) || is.null(init_pi) || is.null(init_nu)) {
      stop("Please provide all initial values of 'init_Sigma', 'init_pi' and 'init_nu'!")
    }
    # Manual initialization provided by user
    if (is.null(init_pi)) pi_k <- rep(1 / K, K) else pi_k <- init_pi
    if (is.null(init_nu)) nu_k <- rep(p + 5, K) else nu_k <- init_nu
    Sigma_k <- init_Sigma
  }

  tau <- matrix(0, n, K)
  loglik_trace <- numeric()

  for (iter in seq_len(maxiter)) {
    # --- 2. E-Step ---
    logdens <- matrix(NA, n, K)

    for (k in seq_len(K)) {
      # Vectorize this loop if possible, but loop is okay for readability
      for (i in seq_len(n)) {
        # if(is.null(Sigma_k[[k]])) browser()
        logdens[i, k] <- dWishart(S_list[[i]], nu_k[k], Sigma_k[[k]],
          logarithm = TRUE
        )
      }
      logdens[, k] <- logdens[, k] + log(pi_k[k])
    }

    # Log-Sum-Exp for numerical stability
    maxlog <- apply(logdens, 1, max)
    # Prevent underflow for very small probabilities
    denom <- maxlog + log(rowSums(exp(logdens - maxlog)))
    logtau <- sweep(logdens, 1, denom, FUN = "-")
    tau <- exp(logtau)

    current_loglik <- sum(denom)
    loglik_trace <- c(loglik_trace, current_loglik)

    if (verbose && iter %% 10 == 0) {
      cat(sprintf(
        "Iter %3d | Loglik: %.4f | Nu: %s\n",
        iter, current_loglik, paste(round(nu_k, 2), collapse = ", ")
      ))
    }

    # Convergence check
    if (iter > 2 && abs(loglik_trace[iter] - loglik_trace[iter - 1]) < tol) {
      if (verbose) cat("Converged at iteration", iter, "\n")
      break
    }

    # --- 3. M-Step ---
    N_k <- colSums(tau)
    pi_k <- N_k / n

    # Prevent collapse of clusters
    if (any(N_k < 1e-6)) {
      # Re-initialize empty cluster to a random data point
      bad_k <- which(N_k < 1e-6)
      for (bk in bad_k) {
        ridx <- sample(n, 1)
        Sigma_k[[bk]] <- S_list[[ridx]] / nu_k[bk]
        N_k[bk] <- 1 # Soft reset
      }
    }

    for (k in seq_len(K)) {
      # Weighted sum of S_i
      # Using Reduce is okay, but loop accumulation is often clearer/faster in R for lists
      Ssum <- matrix(0, p, p)
      for (i in seq_len(n)) {
        if (tau[i, k] > 1e-10) { # Sparse update optimization
          Ssum <- Ssum + tau[i, k] * S_list[[i]]
        }
      }

      # 3a. Update Sigma (conditional on current nu)
      Sigma_k[[k]] <- Ssum / (N_k[k] * nu_k[k])
    }

    for (k in seq_len(K)) {
      # 3b. Update Nu (if requested)
      if (estimate_nu) {
        # Profile Log-Likelihood maximization for nu
        # Target: log(nu/2) + psi(nu/2 + ...) terms vs Data Stats

        T1k <- sum(tau[, k] * logdetS)
        logdetSig <- as.numeric(determinant(Sigma_k[[k]] + diag(p) * 1e-12,
          logarithm = TRUE
        )$modulus)

        lhs <- (T1k / N_k[k]) - logdetSig - p * log(2)

        # Newton-Raphson functions

        f_optim <- function(v) {
          # We want to find root of: LHS - sum(digamma(v/2 + ...))
          # , strictly: E[log|W|] = log|Sigma| + p*log(2) + sum(digamma(...))
          # So sum(digamma(...)) = E[log|W|] - log|Sigma| - p*log(2) = LHS.
          # So we want: LHS - sum(digamma(...)) = 0
          lhs - sum(digamma(v / 2 + (1 - seq_len(p)) / 2))
        }

        fp_optim <- function(v) {
          -0.5 * sum(trigamma(v / 2 + (1 - seq_len(p)) / 2))
        }

        # Newton Iteration
        curr_nu <- nu_k[k]
        for (nm in 1:50) {
          val <- f_optim(curr_nu)
          grad <- fp_optim(curr_nu)
          if (is.na(val) || is.na(grad) || abs(grad) < 1e-8) break

          step <- val / grad
          curr_nu <- curr_nu - step

          # Constraints
          if (curr_nu <= p) curr_nu <- p + 0.1
          if (curr_nu > 1e3) curr_nu <- 1e3
          if (abs(step) < 1e-4) break
        }
        nu_k[k] <- curr_nu

        # 3c. Consistency update: Re-calc Sigma with new nu
        # Sigma_k[[k]] <- Ssum / (N_k[k] * nu_k[k])
      }
    }
    # if (any(nu_k > 1e3)) break
  }

  list(
    pi = pi_k, Sigma = Sigma_k, nu = nu_k, tau = tau,
    loglik = loglik_trace, iterations = iter
  )
}
