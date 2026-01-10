##==========================================================
## This chunck of code is inserted in file 'sim_Figure1A.R'
##==========================================================

pi.error2 <- matrix(0, nrow = n_sim, ncol = 6)
nu.error2 <- matrix(0, nrow = n_sim, ncol = 6)
Sigma.error2 <- matrix(0, nrow = n_sim, ncol = 9)
beta.error2 <- matrix(0, nrow = n_sim, ncol = 6)
colnames(pi.error2) <- colnames(nu.error2) <- colnames(beta.error2) <- c("BayesMean1", "BayesMean2", "BMA1", "BMA2", "EM1", "EM2")
colnames(Sigma.error2) <- c("BayesMean1", "BayesMean2", "BMA1", "BMA2", "EM1", "EM2", "BayesMean1-KL", "BayesMean2-KL", "EM-KL")

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

  # TrueOrder <- c(order(nu.true[-K]), K) # MoE model set last cluster as reference
  # BayesOrder <-  c(order(colMeans(fit$MoEfitBayes$nu)[-K]), K) # MoE model set last cluster as reference
  TrueOrder <- order(nu.true)
  BayesOrder <-  order(colMeans(fit$MoEfitBayes$nu[-c(1:burnin),]))
  pis <- t(apply(fit$MoEfitBayes$z_samples, 1, function(xx) c(mean(xx==1), mean(xx==2), mean(xx==3))))
  pi.error2[i, 1] <- mean(abs(colMeans(pis[-c(1:burnin),])[BayesOrder] - pi.true[TrueOrder]))
  # pi.error2[i, 2] <- sqrt(mean((colMeans(fit$MoEfitBayes$pi)[BayesOrder] - pi.true[TrueOrder])^2))
  pi.error2[i, 2] <- mean((colMeans(pis[-c(1:burnin),])[BayesOrder] - pi.true[TrueOrder])^2)
  # pi.error2[i, 3] <- mean( apply(fit$MoEfitBayes$pi[,BayesOrder], 1, function(xx) mean(abs(xx - pi.true[TrueOrder]))) )
  # # pi.error2[i, 4] <- mean( apply(fit$MoEfitBayes$pi[,BayesOrder], 1, function(xx) sqrt(mean((xx - pi.true[TrueOrder])^2))) )
  # pi.error2[i, 4] <- mean( apply(fit$MoEfitBayes$pi[,BayesOrder], 1, function(xx) mean((xx - pi.true[TrueOrder])^2)) )

  ## the clusters' order by Bayesian model averaging is determined at every MCMC iteration
  for(j in (burnin+1):NROW(fit$MoEfitBayes$nu)) {
    BayesOrderBMA <-  order(fit$MoEfitBayes$nu[j,])
    pi.error2[i, 3] <- pi.error2[i, 3] + mean(abs(pis[j,BayesOrderBMA] - pi.true[TrueOrder]))
    pi.error2[i, 4] <- pi.error2[i, 4] + mean((pis[j,BayesOrderBMA] - pi.true[TrueOrder])^2)

    nu.error2[i, 3] <- nu.error2[i, 3] + mean(abs(fit$MoEfitBayes$nu[j,BayesOrderBMA] - nu.true[TrueOrder]))
    nu.error2[i, 4] <- nu.error2[i, 4] + mean((fit$MoEfitBayes$nu[j,BayesOrderBMA] - nu.true[TrueOrder])^2)

    for(k in 1:K) {
      s0 <- Sigma.true[TrueOrder][[k]]
      s <- fit$MoEfitBayes$Sigma[[j]][,,BayesOrderBMA][,,k]
      tmp <- s[lower.tri(s, diag = TRUE)] - s0[lower.tri(s, diag = TRUE)]
      Sigma.error2[i, 3] <- Sigma.error2[i, 3] + sum(abs(tmp))
      Sigma.error2[i, 4] <- Sigma.error2[i, 4] + sum(tmp^2)

      ## KL divergence
      chol_Sig <- chol(s)
      s_inv <- chol2inv(chol_Sig)
      Sigma.error2[i, 8] <- Sigma.error2[i, 8] + 0.5*( sum(diag(s_inv%*%s0)) - as.numeric(determinant(s_inv%*%s0, logarithm = TRUE)$modulus) - ncol(s0))
    }

    beta_single <- matrix(fit$MoEfitBayes$Beta_samples[j,,], ncol = K)#[, -K] ## not subtract (intercepts and) reference cluster
    if(NROW(beta_single) > 1) {
      # tmp <- beta_single[-1, ] - beta.true[,TrueOrder]#[, -K] #remove intercept from 'beta_single'
      tmp <- beta_bias_func(beta_single[-1, ], beta.true, BayesOrderBMA, TrueOrder)
      beta.error2[i, 3] <- beta.error2[i, 3] + sum(abs(tmp))
      beta.error2[i, 4] <- beta.error2[i, 4] + sum(tmp^2)
    }
  }
  pi.error2[i, 3] <- pi.error2[i, 3] / (NROW(fit$MoEfitBayes$nu) - burnin)
  pi.error2[i, 4] <- pi.error2[i, 4] / (NROW(fit$MoEfitBayes$nu) - burnin)
  nu.error2[i, 3] <- nu.error2[i, 3] / (NROW(fit$MoEfitBayes$nu) - burnin)
  nu.error2[i, 4] <- nu.error2[i, 4] / (NROW(fit$MoEfitBayes$nu) - burnin)
  Sigma.error2[i, 3] <- Sigma.error2[i, 3] / (NROW(fit$MoEfitBayes$nu) - burnin)
  Sigma.error2[i, 4] <- Sigma.error2[i, 4] / (NROW(fit$MoEfitBayes$nu) - burnin)
  Sigma.error2[i, 8] <- Sigma.error2[i, 8] / (NROW(fit$MoEfitBayes$nu) - burnin)
  beta.error2[i, 3] <- beta.error2[i, 3] / (NROW(fit$MoEfitBayes$nu) - burnin)
  beta.error2[i, 4] <- beta.error2[i, 4] / (NROW(fit$MoEfitBayes$nu) - burnin)

  nu.error2[i, 1] <- mean(abs(colMeans(fit$MoEfitBayes$nu[-c(1:burnin),])[BayesOrder] - nu.true[TrueOrder]))
  # nu.error2[i, 2] <- sqrt(mean((colMeans(fit$MoEfitBayes$nu)[BayesOrder] - nu.true[TrueOrder])^2))
  nu.error2[i, 2] <- mean((colMeans(fit$MoEfitBayes$nu[-c(1:burnin),])[BayesOrder] - nu.true[TrueOrder])^2)
  # nu.error2[i, 3] <- mean( apply(fit$MoEfitBayes$nu[,BayesOrder], 1, function(xx) mean(abs(xx - nu.true[TrueOrder]))) )
  # # nu.error2[i, 4] <- mean( apply(fit$MoEfitBayes$nu[,BayesOrder], 1, function(xx) sqrt(mean((xx - nu.true[TrueOrder])^2))) )
  # nu.error2[i, 4] <- mean( apply(fit$MoEfitBayes$nu[,BayesOrder], 1, function(xx) mean((xx - nu.true[TrueOrder])^2)) )

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
    # sqrt( sum((s[lower.tri(s, diag = TRUE)] - s0[lower.tri(s, diag = TRUE)])^2) ) /
    #   sqrt( sum((s0[lower.tri(s, diag = TRUE)])^2) )
  }) )
  # Sigma.error2[i, 3] <- mean( unlist(lapply(fit$MoEfitBayes$Sigma, function(Sigma_single){
  #   mean( sapply(1:K, function(xx) {
  #   s0 <- Sigma.true[TrueOrder][[xx]]
  #   s <- Sigma_single[,,BayesOrder][,,xx]
  #   sum(abs(s[lower.tri(s, diag = TRUE)] - s0[lower.tri(s, diag = TRUE)]))
  #   } ))
  # })) )
  # Sigma.error2[i, 4] <- mean( unlist(lapply(fit$MoEfitBayes$Sigma, function(Sigma_single){
  #   sqrt(mean( sapply(1:K, function(xx) {
  #   s0 <- Sigma.true[TrueOrder][[xx]]
  #   s <- Sigma_single[,,BayesOrder][,,xx]
  #   sum((s[lower.tri(s, diag = TRUE)] - s0[lower.tri(s, diag = TRUE)])^2)
  #   } ) ))
  # })) )
  ## KL divergence
  Sigma.error2[i, 7] <- mean( sapply(1:K, function(xx) {
    s0 <- Sigma.true[TrueOrder][[xx]]
    s <- Sigma.bayes[,,BayesOrder][,,xx]
    chol_Sig <- chol(s)
    s_inv <- chol2inv(chol_Sig)
    0.5*( sum(diag(s_inv%*%s0)) - as.numeric(determinant(s_inv%*%s0, logarithm = TRUE)$modulus) - ncol(s0))
  }) )

  beta.bayes <- apply(fit$MoEfitBayes$Beta_samples, c(2, 3), function(xx){mean(xx[-c(1:burnin)])})[-1, ]#-K] ## subtract intercepts #and reference cluster
  # beta.error2[i, 1] <- mean( sapply(1:K, function(xx) sum(abs(beta.bayes[,BayesOrder][,xx] - beta.true[,TrueOrder][,xx])) ))
  # beta.error2[i, 2] <- mean( sapply(1:K, function(xx) sum((beta.bayes[,BayesOrder][,xx] - beta.true[,TrueOrder][,xx])^2)) )
  if(NROW(beta.bayes) > 1 && !is.null(beta.true)) {
    beta.error2[i, 1] <- mean( sum(abs(beta_bias_func(beta.bayes, beta.true, BayesOrder, TrueOrder))) )
    beta.error2[i, 2] <- mean( sum((beta_bias_func(beta.bayes, beta.true, BayesOrder, TrueOrder))^2) )
  }

  # beta.error2[i, 3] <- mean( apply(fit$MoEfitBayes$Beta_samples, 1, function(xx){
  #   beta_single <- xx[-1, ]
  #   mean( sapply(1:K, function(xx) sum(abs(beta_single[,BayesOrder][,xx] - beta.true[,TrueOrder][,xx])) ))
  # }) )
  # beta.error2[i, 4] <- mean( apply(fit$MoEfitBayes$Beta_samples, 1, function(xx){
  #   beta_single <- xx[-1, ]
  #   mean( sapply(1:K, function(xx) sum((beta_single[,BayesOrder][,xx] - beta.true[,TrueOrder][,xx])^2) ))
  # }) )

  # EMorder <- c(order(fit$MoEfitEM$nu[-K]), K)
  EMorder <- order(fit$MoEfitEM$nu)
  pi.error2[i, 5] <- mean(abs(colMeans(fit$MoEfitEM$gamma)[EMorder] - pi.true[TrueOrder]))
  # pi.error2[i, 6] <- sqrt(mean((fit$MoEfitEM$pi[EMorder] - pi.true[TrueOrder])^2))
  pi.error2[i, 6] <- mean((colMeans(fit$MoEfitEM$gamma)[EMorder] - pi.true[TrueOrder])^2)
  nu.error2[i, 5] <- mean(abs(fit$MoEfitEM$nu[EMorder] - nu.true[TrueOrder]))
  # nu.error2[i, 6] <- sqrt(mean((fit$MoEfitEM$nu[EMorder] - nu.true[TrueOrder])^2))
  nu.error2[i, 6] <- mean((fit$MoEfitEM$nu[EMorder] - nu.true[TrueOrder])^2)
  Sigma.error2[i, 5] <- mean( sapply(1:K, function(xx) {
    s0 <- Sigma.true[TrueOrder][[xx]]
    s <- fit$MoEfitEM$Sigma[EMorder][[xx]]
    sum(abs(s[lower.tri(s, diag = TRUE)] - s0[lower.tri(s, diag = TRUE)])) #/ sum(abs(s0[lower.tri(s, diag = TRUE)]))
  }) )
  Sigma.error2[i, 6] <- mean( sapply(1:K, function(xx) {
    s0 <- Sigma.true[TrueOrder][[xx]]
    s <- fit$MoEfitEM$Sigma[EMorder][[xx]]
    sum((s[lower.tri(s, diag = TRUE)] - s0[lower.tri(s, diag = TRUE)])^2)
    # sqrt( sum((s[lower.tri(s, diag = TRUE)] - s0[lower.tri(s, diag = TRUE)])^2) ) /
    #   sqrt( sum((s0[lower.tri(s, diag = TRUE)])^2) )
  }) )
  ## KL divergence
  Sigma.error2[i, 9] <- mean( sapply(1:K, function(xx) {
    s0 <- Sigma.true[TrueOrder][[xx]]
    s <- fit$MoEfitEM$Sigma[EMorder][[xx]]
    chol_Sig <- chol(s)
    s_inv <- chol2inv(chol_Sig)
    0.5*( sum(diag(s_inv%*%s0)) - as.numeric(determinant(s_inv%*%s0, logarithm = TRUE)$modulus) - ncol(s0))
  }) )

  beta.EM <- fit$MoEfitEM$Beta[-1, ]
  if(NROW(beta.EM) > 1 && !is.null(beta.true)) {
    # beta.error2[i, 5] <- mean( sapply(1:K, function(xx) sum(abs(beta.EM[,EMorder][,xx] - beta.true[,TrueOrder][,xx])) ))
    # beta.error2[i, 6] <- mean( sapply(1:K, function(xx) sum((beta.EM[,EMorder][,xx] - beta.true[,TrueOrder][,xx])^2)) )
    beta.error2[i, 5] <- mean( sum(abs(beta_bias_func(beta.EM, beta.true, EMorder, TrueOrder))) )
    beta.error2[i, 6] <- mean( sum((beta_bias_func(beta.EM, beta.true, EMorder, TrueOrder))^2) )
  }
}

# summary of means
error.mean2 <- 
rbind(
  rbind(apply(pi.error2, 2, mean)[c(1,3,5)], apply(pi.error2, 2, mean)[c(2,4,6)]), 
  rbind(apply(nu.error2, 2, mean)[c(1,3,5)], apply(nu.error2, 2, mean)[c(2,4,6)]),
  rbind(apply(Sigma.error2, 2, mean)[c(1,3,5)], apply(Sigma.error2, 2, mean)[c(2,4,6)]), 
  rbind(apply(beta.error2, 2, mean)[c(1,3,5)], apply(beta.error2, 2, mean)[c(2,4,6)])
)
rownames(error.mean2) <- c("pi-norm1", "pi-norm2", "nu-norm1", "nu-norm2", "Sigma-norm1", "Sigma-norm2", "beta-norm1", "beta-norm2")
error.mean2 <- round(error.mean2, digit = 3)

# summary of sd
error.sd2 <- 
rbind(
  rbind(apply(pi.error2, 2, sd)[c(1,3,5)], apply(pi.error2, 2, sd)[c(2,4,6)]), 
  rbind(apply(nu.error2, 2, sd)[c(1,3,5)], apply(nu.error2, 2, sd)[c(2,4,6)]),
  rbind(apply(Sigma.error2, 2, sd)[c(1,3,5)], apply(Sigma.error2, 2, sd)[c(2,4,6)]),
  rbind(apply(beta.error2, 2, sd)[c(1,3,5)], apply(beta.error2, 2, sd)[c(2,4,6)]) 
)
error.sd2 <- round(error.sd2, digit = 4)

