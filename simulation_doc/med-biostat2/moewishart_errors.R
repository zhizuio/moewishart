# pi errors
pi.error <- matrix(nrow = n_sim, ncol = 6)
nu.error <- matrix(nrow = n_sim, ncol = 6)
Sigma.error <- matrix(nrow = n_sim, ncol = 6)
colnames(pi.error) <- colnames(nu.error) <- colnames(Sigma.error) <- c("BayesMean1", "BayesMean2", "BMA1", "BMA2", "EM1", "EM2")

for(i in 1:n_sim) {
  load(paste0(filename_save, i, ".RData"))
  # cat("sim", i, "Bayesian pi=", colMeans(fit$fitBayes$pi), "; EM pi=", fit$fitEM$pi, "\n")
  # cat("sim", i, "Bayesian nu=", colMeans(fit$fitBayes$nu), "; EM nu=", fit$fitEM$nu, "\n")

  TrueOrder <- order(nu.true)
  BayesOrder <-  order(colMeans(fit$fitBayes$nu))
  pi.error[i, 1] <- mean(abs(colMeans(fit$fitBayes$pi)[BayesOrder] - pi.true[TrueOrder]))
  # pi.error[i, 2] <- sqrt(mean((colMeans(fit$fitBayes$pi)[BayesOrder] - pi.true[TrueOrder])^2))
  pi.error[i, 2] <- mean((colMeans(fit$fitBayes$pi)[BayesOrder] - pi.true[TrueOrder])^2)
  pi.error[i, 3] <- mean( apply(fit$fitBayes$pi[,BayesOrder], 1, function(xx) mean(abs(xx - pi.true[TrueOrder]))) )
  # pi.error[i, 4] <- mean( apply(fit$fitBayes$pi[,BayesOrder], 1, function(xx) sqrt(mean((xx - pi.true[TrueOrder])^2))) )
  pi.error[i, 4] <- mean( apply(fit$fitBayes$pi[,BayesOrder], 1, function(xx) mean((xx - pi.true[TrueOrder])^2)) )

  nu.error[i, 1] <- mean(abs(colMeans(fit$fitBayes$nu)[BayesOrder] - nu.true[TrueOrder]))
  # nu.error[i, 2] <- sqrt(mean((colMeans(fit$fitBayes$nu)[BayesOrder] - nu.true[TrueOrder])^2))
  nu.error[i, 2] <- mean((colMeans(fit$fitBayes$nu)[BayesOrder] - nu.true[TrueOrder])^2)
  nu.error[i, 3] <- mean( apply(fit$fitBayes$nu[,BayesOrder], 1, function(xx) mean(abs(xx - nu.true[TrueOrder]))) )
  # nu.error[i, 4] <- mean( apply(fit$fitBayes$nu[,BayesOrder], 1, function(xx) sqrt(mean((xx - nu.true[TrueOrder])^2))) )
  nu.error[i, 4] <- mean( apply(fit$fitBayes$nu[,BayesOrder], 1, function(xx) mean((xx - nu.true[TrueOrder])^2)) )

  if (is.list(fit$fitBayes$Sigma)) {
    fit$fitBayes$Sigma <- lapply(fit$fitBayes$Sigma, function(xx) simplify2array(xx))
  }
  Sigma.bayes <- Reduce("+", fit$fitBayes$Sigma)/length(fit$fitBayes$Sigma)
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
  Sigma.error[i, 3] <- mean( unlist(lapply(fit$fitBayes$Sigma, function(Sigma_single){
    mean( sapply(1:K, function(xx) {
    s0 <- Sigma.true[TrueOrder][[xx]]
    s <- Sigma_single[,,BayesOrder][,,xx]
    sum(abs(s[lower.tri(s, diag = TRUE)] - s0[lower.tri(s, diag = TRUE)]))
    } ))
  })) )
  # Sigma.error[i, 4] <- mean( unlist(lapply(fit$fitBayes$Sigma, function(Sigma_single){
  #   sqrt(mean( sapply(1:K, function(xx) (Sigma_single[,,BayesOrder][,,xx] - Sigma.true[TrueOrder][[xx]])^2) ))
  # })) )
  Sigma.error[i, 4] <- mean( unlist(lapply(fit$fitBayes$Sigma, function(Sigma_single){
    sqrt(mean( sapply(1:K, function(xx) {
    s0 <- Sigma.true[TrueOrder][[xx]]
    s <- Sigma_single[,,BayesOrder][,,xx]
    sum((s[lower.tri(s, diag = TRUE)] - s0[lower.tri(s, diag = TRUE)])^2)
    } ) ))
  })) )

  EMorder <- order(fit$fitEM$nu)
  pi.error[i, 5] <- mean(abs(fit$fitEM$pi[EMorder] - pi.true[TrueOrder]))
  # pi.error[i, 6] <- sqrt(mean((fit$fitEM$pi[EMorder] - pi.true[TrueOrder])^2))
  pi.error[i, 6] <- mean((fit$fitEM$pi[EMorder] - pi.true[TrueOrder])^2)
  nu.error[i, 5] <- mean(abs(fit$fitEM$nu[EMorder] - nu.true[TrueOrder]))
  # nu.error[i, 6] <- sqrt(mean((fit$fitEM$nu[EMorder] - nu.true[TrueOrder])^2))
  nu.error[i, 6] <- mean((fit$fitEM$nu[EMorder] - nu.true[TrueOrder])^2)
  Sigma.error[i, 5] <- mean( sapply(1:K, function(xx) {
    s0 <- Sigma.true[TrueOrder][[xx]]
    s <- fit$fitEM$Sigma[EMorder][[xx]]
    sum(abs(s[lower.tri(s, diag = TRUE)] - s0[lower.tri(s, diag = TRUE)]))
  }) )
  # Sigma.error[i, 6] <- sqrt(mean( sapply(1:K, function(xx) sum((fit$fitEM$Sigma[EMorder][[xx]] - Sigma.true[TrueOrder][[xx]])^2)) ))
  Sigma.error[i, 6] <- mean( sapply(1:K, function(xx) {
    s0 <- Sigma.true[TrueOrder][[xx]]
    s <- fit$fitEM$Sigma[EMorder][[xx]]
    sum((s[lower.tri(s, diag = TRUE)] - s0[lower.tri(s, diag = TRUE)])^2)
  }) )
}

# summary of means
error.mean <- 
rbind(
  rbind(apply(pi.error, 2, mean)[c(1,3,5)], apply(pi.error, 2, mean)[c(2,4,6)]), 
  rbind(apply(nu.error, 2, mean)[c(1,3,5)], apply(nu.error, 2, mean)[c(2,4,6)]),
  rbind(apply(Sigma.error, 2, mean)[c(1,3,5)], apply(Sigma.error, 2, mean)[c(2,4,6)]) 
)
rownames(error.mean) <- c("pi-norm1", "pi-norm2", "nu-norm1", "nu-norm2", "Sigma-norm1", "Sigma-norm2")
error.mean <- round(error.mean, digit = 3)

# summary of sd
error.sd <- 
rbind(
  rbind(apply(pi.error, 2, sd)[c(1,3,5)], apply(pi.error, 2, sd)[c(2,4,6)]), 
  rbind(apply(nu.error, 2, sd)[c(1,3,5)], apply(nu.error, 2, sd)[c(2,4,6)]),
  rbind(apply(Sigma.error, 2, sd)[c(1,3,5)], apply(Sigma.error, 2, sd)[c(2,4,6)]) 
)
error.sd <- round(error.sd, digit = 4)

