##==========================================================
## This chunck of code is inserted in file 'sim_Figure1A.R'
##==========================================================

pi.error2 <- matrix(0, nrow = n_sim, ncol = 4)
nu.error2 <- matrix(0, nrow = n_sim, ncol = 4)
Sigma.error2 <- matrix(0, nrow = n_sim, ncol = 4)
beta.error2 <- matrix(0, nrow = n_sim, ncol = 4)
colnames(pi.error2) <- colnames(nu.error2) <- colnames(beta.error2) <- colnames(Sigma.error2) <- c("BayesMean1", "BayesMean2", "EM1", "EM2")

## Function to compute the bias of coefficients with matched reference cluster
beta_bias_func <- function(B_hat, B_true, ord_est, ord_true) {
  # B_hat: p x K estimated coefficients with baseline column = 0 at column K (before matching)
  # perm: permutation such that estimated column perm[j] matches true column j
  # Reorder estimated columns into true order:

  # perm[j] = estimated index matched to true index j
  q <- NCOL(B_true)
  perm <- numeric(q)
  perm[ord_true] <- ord_est
  # perm <- as.integer(perm)

  # Sanity check: after reordering, the columns of B_hat should correspond to nu_true order
  # nu_est[perm] should be aligned in rank with nu_true
  # print(nu_true); print(nu_est[perm])

  # 2) Reorder estimated coefficients to the true ordering
  B_hat_ord <- B_hat[, perm, drop = FALSE]

  # 3) Rebase so that the matched true baseline column j_true is zero in the estimated matrix
  j_true <- which(colSums(B_true == 0) == q)
  b_base <- B_hat_ord[, j_true, drop = FALSE]   # p x 1
  B_hat_rebased <- sweep(B_hat_ord, 1, b_base, FUN = "-")        # subtract baseline column from all columns

  # 4) Compute bias of coefficients (now baseline-aligned)
  bias_mat <- B_hat_rebased - B_true

  return(bias_mat)
}

for(i in 1:n_sim) {
  load(paste0(filename_save, i, ".RData"))

  TrueOrder <- order(nu.true)
  BayesOrder <-  order(colMeans(fit$MoEfitBayes$nu[-c(1:burnin),]))
  pis <- t(apply(fit$MoEfitBayes$z_samples, 1, function(xx) c(mean(xx==1), mean(xx==2), mean(xx==3))))
  pi.error2[i, 1] <- mean(abs(colMeans(pis[-c(1:burnin),])[BayesOrder] - pi.true[TrueOrder]))
  pi.error2[i, 2] <- mean((colMeans(pis[-c(1:burnin),])[BayesOrder] - pi.true[TrueOrder])^2)

  nu.error2[i, 1] <- mean(abs(colMeans(fit$MoEfitBayes$nu[-c(1:burnin),])[BayesOrder] - nu.true[TrueOrder]))
  nu.error2[i, 2] <- mean((colMeans(fit$MoEfitBayes$nu[-c(1:burnin),])[BayesOrder] - nu.true[TrueOrder])^2)

  if (is.list(fit$MoEfitBayes$Sigma)) {
    fit$MoEfitBayes$Sigma <- lapply(fit$MoEfitBayes$Sigma[-c(1:burnin)], function(xx) simplify2array(xx))
  }
  Sigma.bayes <- Reduce("+", fit$MoEfitBayes$Sigma)/length(fit$MoEfitBayes$Sigma)
  Sigma.error2[i, 1] <- mean( sapply(1:K, function(xx) {
    s0 <- Sigma.true[TrueOrder][[xx]]
    s <- Sigma.bayes[,,BayesOrder][,,xx]
    sum(abs(s[lower.tri(s,diag=TRUE)] - s0[lower.tri(s0,diag=TRUE)])) #/ sum(abs(s0[lower.tri(s0,diag=TRUE)]))
  } ))
  Sigma.error2[i, 2] <- mean( sapply(1:K, function(xx) {
    s0 <- Sigma.true[TrueOrder][[xx]]
    s <- Sigma.bayes[,,BayesOrder][,,xx]
    sum((s[lower.tri(s, diag = TRUE)] - s0[lower.tri(s, diag = TRUE)])^2)
  }) )

  beta.bayes <- apply(fit$MoEfitBayes$Beta_samples, c(2, 3), function(xx){mean(xx[-c(1:burnin)])})[-1, ]
  if(NROW(beta.bayes) > 1 && !is.null(beta.true)) {
    beta.error2[i, 1] <- mean( sum(abs(beta_bias_func(beta.bayes, beta.true, BayesOrder, TrueOrder))) )
    beta.error2[i, 2] <- mean( sum((beta_bias_func(beta.bayes, beta.true, BayesOrder, TrueOrder))^2) )
  }

  EMorder <- order(fit$MoEfitEM$nu)
  pi.error2[i, 3] <- mean(abs(colMeans(fit$MoEfitEM$gamma)[EMorder] - pi.true[TrueOrder]))
  pi.error2[i, 4] <- mean((colMeans(fit$MoEfitEM$gamma)[EMorder] - pi.true[TrueOrder])^2)
  nu.error2[i, 3] <- mean(abs(fit$MoEfitEM$nu[EMorder] - nu.true[TrueOrder]))
  nu.error2[i, 4] <- mean((fit$MoEfitEM$nu[EMorder] - nu.true[TrueOrder])^2)
  Sigma.error2[i, 3] <- mean( sapply(1:K, function(xx) {
    s0 <- Sigma.true[TrueOrder][[xx]]
    s <- fit$MoEfitEM$Sigma[EMorder][[xx]]
    sum(abs(s[lower.tri(s, diag = TRUE)] - s0[lower.tri(s, diag = TRUE)])) 
  }) )
  Sigma.error2[i, 4] <- mean( sapply(1:K, function(xx) {
    s0 <- Sigma.true[TrueOrder][[xx]]
    s <- fit$MoEfitEM$Sigma[EMorder][[xx]]
    sum((s[lower.tri(s, diag = TRUE)] - s0[lower.tri(s, diag = TRUE)])^2)
  }) )

  beta.EM <- fit$MoEfitEM$Beta[-1, ]
  if(NROW(beta.EM) > 1 && !is.null(beta.true)) {
    beta.error2[i, 3] <- mean( sum(abs(beta_bias_func(beta.EM, beta.true, EMorder, TrueOrder))) )
    beta.error2[i, 4] <- mean( sum((beta_bias_func(beta.EM, beta.true, EMorder, TrueOrder))^2) )
  }
}

# summary of means
error.mean2 <- 
rbind(
  rbind(apply(pi.error2, 2, mean)[c(1,3)], apply(pi.error2, 2, mean)[c(2,4)]), 
  rbind(apply(nu.error2, 2, mean)[c(1,3)], apply(nu.error2, 2, mean)[c(2,4)]),
  rbind(apply(Sigma.error2, 2, mean)[c(1,3)], apply(Sigma.error2, 2, mean)[c(2,4)]), 
  rbind(apply(beta.error2, 2, mean)[c(1,3)], apply(beta.error2, 2, mean)[c(2,4)])
)
rownames(error.mean2) <- c("pi-norm1", "pi-norm2", "nu-norm1", "nu-norm2", "Sigma-norm1", "Sigma-norm2", "beta-norm1", "beta-norm2")
error.mean2 <- round(error.mean2, digit = 3)

# summary of sd
error.sd2 <- 
rbind(
  rbind(apply(pi.error2, 2, sd)[c(1,3)], apply(pi.error2, 2, sd)[c(2,4)]), 
  rbind(apply(nu.error2, 2, sd)[c(1,3)], apply(nu.error2, 2, sd)[c(2,4)]),
  rbind(apply(Sigma.error2, 2, sd)[c(1,3)], apply(Sigma.error2, 2, sd)[c(2,4)]),
  rbind(apply(beta.error2, 2, sd)[c(1,3)], apply(beta.error2, 2, sd)[c(2,4)]) 
)
error.sd2 <- round(error.sd2, digit = 4)

