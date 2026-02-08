##==========================================================
## This chunck of code is inserted in file 'sim_Figure1A.R'
##==========================================================

pi.error <- matrix(0, nrow = n_sim, ncol = 4)
nu.error <- matrix(0, nrow = n_sim, ncol = 4)
Sigma.error <- matrix(0, nrow = n_sim, ncol = 4)
colnames(pi.error) <- colnames(nu.error) <- colnames(Sigma.error) <- c("BayesMean1", "BayesMean2", "EM1", "EM2")

for(i in 1:n_sim) {
  load(paste0(filename_save, i, ".RData"))

  TrueOrder <- order(nu.true)
  BayesOrder <-  order(colMeans(fit$fitBayes$nu[-c(1:burnin),]))
  pi.error[i, 1] <- mean(abs(colMeans(fit$fitBayes$pi[-c(1:burnin),])[BayesOrder] - pi.true[TrueOrder]))
  pi.error[i, 2] <- mean((colMeans(fit$fitBayes$pi[-c(1:burnin),])[BayesOrder] - pi.true[TrueOrder])^2)

  nu.error[i, 1] <- mean(abs(colMeans(fit$fitBayes$nu[-c(1:burnin),])[BayesOrder] - nu.true[TrueOrder]))
  nu.error[i, 2] <- mean((colMeans(fit$fitBayes$nu[-c(1:burnin),])[BayesOrder] - nu.true[TrueOrder])^2)

  if (is.list(fit$fitBayes$Sigma)) {
    fit$fitBayes$Sigma <- lapply(fit$fitBayes$Sigma[-c(1:burnin)], function(xx) simplify2array(xx))
  }
  Sigma.bayes <- Reduce("+", fit$fitBayes$Sigma)/length(fit$fitBayes$Sigma)
  Sigma.error[i, 1] <- mean( sapply(1:K, function(xx) {
    s0 <- Sigma.true[TrueOrder][[xx]]
    s <- Sigma.bayes[,,BayesOrder][,,xx]
    sum(abs(s[lower.tri(s,diag=TRUE)] - s0[lower.tri(s0,diag=TRUE)])) #/ sum(abs(s0[lower.tri(s0,diag=TRUE)]))
  } ))
  Sigma.error[i, 2] <- mean( sapply(1:K, function(xx) {
    s0 <- Sigma.true[TrueOrder][[xx]]
    s <- Sigma.bayes[,,BayesOrder][,,xx]
    sum((s[lower.tri(s, diag = TRUE)] - s0[lower.tri(s, diag = TRUE)])^2)
  }) )

  EMorder <- order(fit$fitEM$nu)
  pi.error[i, 3] <- mean(abs(fit$fitEM$pi[EMorder] - pi.true[TrueOrder]))
  pi.error[i, 4] <- mean((fit$fitEM$pi[EMorder] - pi.true[TrueOrder])^2)
  nu.error[i, 3] <- mean(abs(fit$fitEM$nu[EMorder] - nu.true[TrueOrder]))
  nu.error[i, 4] <- mean((fit$fitEM$nu[EMorder] - nu.true[TrueOrder])^2)
  Sigma.error[i, 3] <- mean( sapply(1:K, function(xx) {
    s0 <- Sigma.true[TrueOrder][[xx]]
    s <- fit$fitEM$Sigma[EMorder][[xx]]
    sum(abs(s[lower.tri(s, diag = TRUE)] - s0[lower.tri(s, diag = TRUE)])) 
  }) )
  Sigma.error[i, 4] <- mean( sapply(1:K, function(xx) {
    s0 <- Sigma.true[TrueOrder][[xx]]
    s <- fit$fitEM$Sigma[EMorder][[xx]]
    sum((s[lower.tri(s, diag = TRUE)] - s0[lower.tri(s, diag = TRUE)])^2)
  }) )
}

# summary of means
error.mean <- 
rbind(
  rbind(apply(pi.error, 2, mean)[c(1,3)], apply(pi.error, 2, mean)[c(2,4)]), 
  rbind(apply(nu.error, 2, mean)[c(1,3)], apply(nu.error, 2, mean)[c(2,4)]),
  rbind(apply(Sigma.error, 2, mean)[c(1,3)], apply(Sigma.error, 2, mean)[c(2,4)]) 
)
rownames(error.mean) <- c("pi-norm1", "pi-norm2", "nu-norm1", "nu-norm2", "Sigma-norm1", "Sigma-norm2")
error.mean <- round(error.mean, digit = 3)

# summary of sd
error.sd <- 
rbind(
  rbind(apply(pi.error, 2, sd)[c(1,3)], apply(pi.error, 2, sd)[c(2,4)]), 
  rbind(apply(nu.error, 2, sd)[c(1,3)], apply(nu.error, 2, sd)[c(2,4)]),
  rbind(apply(Sigma.error, 2, sd)[c(1,3)], apply(Sigma.error, 2, sd)[c(2,4)]) 
)
error.sd <- round(error.sd, digit = 4)

