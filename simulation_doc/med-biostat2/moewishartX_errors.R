
# errors
pi.error <- matrix(nrow = n_sim, ncol = 6)
nu.error <- matrix(nrow = n_sim, ncol = 6)
Sigma.error <- matrix(nrow = n_sim, ncol = 6)
beta.error <- matrix(nrow = n_sim, ncol = 6)
colnames(pi.error) <- colnames(nu.error) <- colnames(Sigma.error) <- colnames(beta.error) <- 
c("BayesMean1", "BayesMean2", "BMA1", "BMA2", "EM1", "EM2")

for(i in 1:n_sim) {
  load(paste0(filename_save, i, ".RData"))

  # TrueOrder <- c(order(nu.true[-K]), K) # MoE model set last cluster as reference
  # BayesOrder <-  c(order(colMeans(fit$MoEfitBayes$nu)[-K]), K) # MoE model set last cluster as reference
  TrueOrder <- order(nu.true)
  BayesOrder <-  order(colMeans(fit$MoEfitBayes$nu))
  pi.error[i, 1] <- mean(abs(colMeans(fit$MoEfitBayes$pi)[BayesOrder] - pi.true[TrueOrder]))
  # pi.error[i, 2] <- sqrt(mean((colMeans(fit$MoEfitBayes$pi)[BayesOrder] - pi.true[TrueOrder])^2))
  pi.error[i, 2] <- mean((colMeans(fit$MoEfitBayes$pi)[BayesOrder] - pi.true[TrueOrder])^2)
  pi.error[i, 3] <- mean( apply(fit$MoEfitBayes$pi[,BayesOrder], 1, function(xx) mean(abs(xx - pi.true[TrueOrder]))) )
  # pi.error[i, 4] <- mean( apply(fit$MoEfitBayes$pi[,BayesOrder], 1, function(xx) sqrt(mean((xx - pi.true[TrueOrder])^2))) )
  pi.error[i, 4] <- mean( apply(fit$MoEfitBayes$pi[,BayesOrder], 1, function(xx) mean((xx - pi.true[TrueOrder])^2)) )

  nu.error[i, 1] <- mean(abs(colMeans(fit$MoEfitBayes$nu)[BayesOrder] - nu.true[TrueOrder]))
  # nu.error[i, 2] <- sqrt(mean((colMeans(fit$MoEfitBayes$nu)[BayesOrder] - nu.true[TrueOrder])^2))
  nu.error[i, 2] <- mean((colMeans(fit$MoEfitBayes$nu)[BayesOrder] - nu.true[TrueOrder])^2)
  nu.error[i, 3] <- mean( apply(fit$MoEfitBayes$nu[,BayesOrder], 1, function(xx) mean(abs(xx - nu.true[TrueOrder]))) )
  # nu.error[i, 4] <- mean( apply(fit$MoEfitBayes$nu[,BayesOrder], 1, function(xx) sqrt(mean((xx - nu.true[TrueOrder])^2))) )
  nu.error[i, 4] <- mean( apply(fit$MoEfitBayes$nu[,BayesOrder], 1, function(xx) mean((xx - nu.true[TrueOrder])^2)) )

  if (is.list(fit$MoEfitBayes$Sigma)) {
    fit$MoEfitBayes$Sigma <- lapply(fit$MoEfitBayes$Sigma, function(xx) simplify2array(xx))
  }
  Sigma.bayes <- Reduce("+", fit$MoEfitBayes$Sigma)/length(fit$MoEfitBayes$Sigma)
  Sigma.error[i, 1] <- mean( sapply(1:K, function(xx) {
    s0 <- Sigma.true[TrueOrder][[xx]]
    s <- Sigma.bayes[,,BayesOrder][,,xx]
    sum(abs(s[lower.tri(s,diag=TRUE)] - s0[lower.tri(s0,diag=TRUE)]))
  } ))
  # Sigma.error[i, 2] <- sqrt(mean( sapply(1:K, function(xx) sum((Sigma.bayes[,,BayesOrder][,,xx] - Sigma.true[TrueOrder][[xx]])^2)) ))
  # Sigma.error[i, 2] <- sqrt(mean( sapply(1:K, function(xx) (Sigma.bayes[,,BayesOrder][,,xx] - Sigma.true[TrueOrder][[xx]])^2) ))
  Sigma.error[i, 2] <- mean( sapply(1:K, function(xx) {
    s0 <- Sigma.true[TrueOrder][[xx]]
    s <- Sigma.bayes[,,BayesOrder][,,xx]
    sum((s[lower.tri(s, diag = TRUE)] - s0[lower.tri(s, diag = TRUE)])^2)
  }) )
  Sigma.error[i, 3] <- mean( unlist(lapply(fit$MoEfitBayes$Sigma, function(Sigma_single){
    mean( sapply(1:K, function(xx) {
    s0 <- Sigma.true[TrueOrder][[xx]]
    s <- Sigma_single[,,BayesOrder][,,xx]
    sum(abs(s[lower.tri(s, diag = TRUE)] - s0[lower.tri(s, diag = TRUE)]))
    } ))
  })) )
  # Sigma.error[i, 4] <- mean( unlist(lapply(fit$MoEfitBayes$Sigma, function(Sigma_single){
  #   sqrt(mean( sapply(1:K, function(xx) (Sigma_single[,,BayesOrder][,,xx] - Sigma.true[TrueOrder][[xx]])^2) ))
  # })) )
  Sigma.error[i, 4] <- mean( unlist(lapply(fit$MoEfitBayes$Sigma, function(Sigma_single){
    sqrt(mean( sapply(1:K, function(xx) {
    s0 <- Sigma.true[TrueOrder][[xx]]
    s <- Sigma_single[,,BayesOrder][,,xx]
    sum((s[lower.tri(s, diag = TRUE)] - s0[lower.tri(s, diag = TRUE)])^2)
    } ) ))
  })) )

  beta.bayes <- apply(fit$MoEfitBayes$Beta_samples, c(2, 3), mean)[-1, ]#-K] ## subtract intercepts #and reference cluster
  beta.error[i, 1] <- mean( sapply(1:K, function(xx) sum(abs(beta.bayes[,BayesOrder][,xx] - beta.true[,TrueOrder][,xx])) ))
  beta.error[i, 2] <- mean( sapply(1:K, function(xx) sum((beta.bayes[,BayesOrder][,xx] - beta.bayes[,TrueOrder][,xx])^2)) )
  beta.error[i, 3] <- mean( apply(fit$MoEfitBayes$Beta_samples, 1, function(xx){
    beta_single <- xx[-1, ]
    mean( sapply(1:K, function(xx) sum(abs(beta_single[,BayesOrder][,xx] - beta.true[,TrueOrder][,xx])) ))
  }) )
  beta.error[i, 4] <- mean( apply(fit$MoEfitBayes$Beta_samples, 1, function(xx){
    beta_single <- xx[-1, ]
    mean( sapply(1:K, function(xx) sum((beta_single[,BayesOrder][,xx] - beta.true[,TrueOrder][,xx])^2) ))
  }) )

  # EMorder <- c(order(fit$MoEfitEM$nu[-K]), K)
  EMorder <- order(fit$MoEfitEM$nu)
  pi.error[i, 5] <- mean(abs(fit$MoEfitEM$gamma[EMorder] - pi.true[TrueOrder]))
  # pi.error[i, 6] <- sqrt(mean((fit$MoEfitEM$pi[EMorder] - pi.true[TrueOrder])^2))
  pi.error[i, 6] <- mean((fit$MoEfitEM$gamma[EMorder] - pi.true[TrueOrder])^2)
  nu.error[i, 5] <- mean(abs(fit$MoEfitEM$nu[EMorder] - nu.true[TrueOrder]))
  # nu.error[i, 6] <- sqrt(mean((fit$MoEfitEM$nu[EMorder] - nu.true[TrueOrder])^2))
  nu.error[i, 6] <- mean((fit$MoEfitEM$nu[EMorder] - nu.true[TrueOrder])^2)
  Sigma.error[i, 5] <- mean( sapply(1:K, function(xx) {
    s0 <- Sigma.true[TrueOrder][[xx]]
    s <- fit$MoEfitEM$Sigma[EMorder][[xx]]
    sum(abs(s[lower.tri(s, diag = TRUE)] - s0[lower.tri(s, diag = TRUE)]))
  }) )
  # Sigma.error[i, 6] <- sqrt(mean( sapply(1:K, function(xx) sum((fit$MoEfitEM$Sigma[EMorder][[xx]] - Sigma.true[TrueOrder][[xx]])^2)) ))
  Sigma.error[i, 6] <- mean( sapply(1:K, function(xx) {
    s0 <- Sigma.true[TrueOrder][[xx]]
    s <- fit$MoEfitEM$Sigma[EMorder][[xx]]
    sum((s[lower.tri(s, diag = TRUE)] - s0[lower.tri(s, diag = TRUE)])^2)
  }) )
  beta.EM <- fit$MoEfitEM$Beta[-1, ]
  beta.error[i, 5] <- mean( sapply(1:K, function(xx) sum(abs(beta.EM[,EMorder][,xx] - beta.true[,TrueOrder][,xx])) ))
  beta.error[i, 6] <- mean( sapply(1:K, function(xx) sum((beta.EM[,EMorder][,xx] - beta.bayes[,TrueOrder][,xx])^2)) )

}

# summary of means
error.mean2 <- 
rbind(
  rbind(apply(pi.error, 2, mean)[c(1,3,5)], apply(pi.error, 2, mean)[c(2,4,6)]), 
  rbind(apply(nu.error, 2, mean)[c(1,3,5)], apply(nu.error, 2, mean)[c(2,4,6)]),
  rbind(apply(Sigma.error, 2, mean)[c(1,3,5)], apply(Sigma.error, 2, mean)[c(2,4,6)]), 
  rbind(apply(beta.error, 2, mean)[c(1,3,5)], apply(beta.error, 2, mean)[c(2,4,6)])
)
rownames(error.mean2) <- c("pi-norm1", "pi-norm2", "nu-norm1", "nu-norm2", "Sigma-norm1", "Sigma-norm2", "beta-norm1", "beta-norm2")
error.mean2 <- round(error.mean2, digit = 3)

# summary of sd
error.sd2 <- 
rbind(
  rbind(apply(pi.error, 2, sd)[c(1,3,5)], apply(pi.error, 2, sd)[c(2,4,6)]), 
  rbind(apply(nu.error, 2, sd)[c(1,3,5)], apply(nu.error, 2, sd)[c(2,4,6)]),
  rbind(apply(Sigma.error, 2, sd)[c(1,3,5)], apply(Sigma.error, 2, sd)[c(2,4,6)]),
  rbind(apply(beta.error, 2, sd)[c(1,3,5)], apply(beta.error, 2, sd)[c(2,4,6)]) 
)
error.sd2 <- round(error.sd2, digit = 4)

