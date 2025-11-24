#' @title Optimized Gibbs Sampler
#'
#' @description
#' TBA
#'
#' @name mewishart
#'
#' @importFrom utils combn
#' @importFrom stats kmeans
#'
#' @param S_list TBA
#' @param K TBA
#' @param niter TBA
#' @param burnin TBA
#' @param thin TBA
#' @param alpha TBA
#' @param nu0 TBA
#' @param Psi0 TBA
#' @param init_nu TBA
#' @param sample_nu TBA
#' @param nu_prior_a TBA
#' @param nu_prior_b TBA
#' @param mh_sigma TBA
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
mewishart <- function(S_list,
                      K,
                      niter = 3000,
                      burnin = 1000,
                      thin = 1,
                      alpha = NULL,
                      nu0 = NULL,
                      Psi0 = NULL,
                      init_nu = NULL,
                      sample_nu = TRUE,
                      nu_prior_a = 2,
                      nu_prior_b = 0.1,
                      mh_sigma = 0.1,
                      verbose = TRUE) {
  # -- 1. Pre-processing and Pre-allocation --
  n <- length(S_list)
  p <- nrow(S_list[[1]])

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

  pi_k <- table(factor(z, levels = 1:K)) / n
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
  nsave <- floor((niter - burnin) / thin)
  if (nsave < 1) nsave <- 1
  out_pi <- matrix(NA, nrow = nsave, ncol = K)
  out_nu <- matrix(NA, nrow = nsave, ncol = K)
  out_Sigma <- vector("list", nsave)
  out_z <- matrix(NA, nrow = nsave, ncol = n)
  logliks <- numeric(niter)
  iter_save <- 0

  # Pre-allocate reusable vectors
  logpost <- matrix(0, n, K)

  # Start Timer
  start_time <- Sys.time()

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

      logpost[, k] <- log(pi_k[k] + 1e-300) + term1 + term2 + term3 + term4
    }

    # Sample z
    # Vectorized sampling is hard in base R, looping sample.int is okay
    # but we can optimize the probability normalization
    # Using pure R loop for sampling is usually fast enough compared to the math above
    for (i in 1:n) {
      lp <- logpost[i, ]
      lp <- lp - max(lp)
      prob <- exp(lp)
      z[i] <- sample.int(K, 1, prob = prob)
    }

    # --- Step 2: Update Weights pi ---
    n_k <- as.numeric(table(factor(z, levels = 1:K)))
    pi_k <- as.numeric(rdirichlet(1, alpha + n_k))

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
        Psi_post <- Psi0 + S_sum

        # Invert Psi_post once for sampling
        # Use tryCatch for numerical stability
        Psi_post_inv <- tryCatch(solve(Psi_post), error = function(e) solve(Psi_post + diag(1e-6, p)))

        Sigma_k[, , k] <- sampleIW(nu_post, Psi_post_inv)
      }
    }

    # --- Step 4: Update nu (MH) ---
    if (sample_nu) {
      for (k in 1:K) {
        curr_nu <- nu_k[k]
        prop_log <- rnorm(1, log(curr_nu), mh_sigma)
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
            Sig_inv <- chol2inv(chol_Sig)

            # Sum of traces = Tr(Sig^-1 * Sum(S))
            # We can calculate Sum(S) fast
            S_sum_vec <- colSums(S_mat[idx, , drop = FALSE])
            sum_tr_val <- sum(as.vector(Sig_inv) * S_sum_vec)

            # Define a mini function for log-lik given Sufficient Stats
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

          # Priors + Jacobian
          lp_old <- (nu_prior_a - 1) * log(curr_nu) - nu_prior_b * curr_nu + log(curr_nu)
          lp_new <- (nu_prior_a - 1) * log(prop_nu) - nu_prior_b * prop_nu + log(prop_nu)

          if (log(runif(1)) < (ll_new + lp_new) - (ll_old + lp_old)) {
            nu_k[k] <- prop_nu
          }
        }
      }
    }

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

    # --- Save ---
    if (iter > burnin && ((iter - burnin) %% thin == 0)) {
      iter_save <- iter_save + 1
      if (iter_save <= nsave) {
        out_pi[iter_save, ] <- pi_k
        out_nu[iter_save, ] <- nu_k
        out_Sigma[[iter_save]] <- Sigma_k
        out_z[iter_save, ] <- z
      }
    }

    if (verbose && (iter %% 500 == 0 || iter == 1)) {
      # Calculate speed
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      rate <- iter / elapsed
      cat(sprintf(
        "Iter %4d | LL=%.1f | %.1f iter/sec\n",
        iter, logliks[iter], rate
      ))
    }
  }

  sigma_posterior_mean <- Reduce("+", out_Sigma) / length(out_Sigma)

  list(
    pi = out_pi, nu = out_nu, Sigma = out_Sigma, z = out_z,
    sigma_posterior_mean = sigma_posterior_mean, loglik = logliks
  )
}
