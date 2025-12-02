rm(list = ls())
filename_save <- "/data/zhiz/moewishart/results/n200_p2K3_sim"
args <- commandArgs(trailingOnly = TRUE)
seed.idx <- as.numeric(args)

con <- file(paste0(filename_save, seed.idx,".log"))
sink(con, append=TRUE)
sink(con, append=TRUE, type="message")


n <- 200 # subjects
p <- 2

library(moewishart)

set.seed(seed.idx)
dat <- simData(n, p,
  pis = c(0.35, 0.40, 0.25),
  nus = c(8, 12, 3),
  Sigma = NULL # default Sigmas pre-defined in function 'simData()'
)
S_list <- dat$S
Sigma_list <- dat$Sigma_list
my_K_intial <- 3
init_pi <- c(0.33, 0.33, 0.34)

# run Bayesian mixture model
set.seed(123)
fitBayes <- moewishart(S_list,
  K = my_K_intial, #init_pi = init_pi, 
  mh_sigma = 0.2,
  niter = 10000, burnin = 1000, thin = 1, verbose = TRUE
)

colMeans(fitBayes$pi) # estimated pi
colMeans(dat$pi) # true pi

apply(fitBayes$nu, 2, mean) # estimated nu
dat$nu # true nu


# run mixture model with EM algorithm
set.seed(123)
fitEM <- moewishart(S_list,
  method = "em",
  K = my_K_intial, #init_pi = init_pi, 
  niter = 10000, verbose = TRUE
)
# EM estimates
fitEM$pi
colMeans(dat$pi) # true pi

fitEM$nu
dat$nu # true nu


# run Bayesian mixture model
set.seed(123)
MoEfitBayes <- moewishartX(
  S_list, X = matrix(rep(1, n), ncol = 1),
  K = my_K_intial, #init_pi = init_pi, 
  mh_sigma = 0.2,
  niter = 10000, burnin = 1000, thin = 1, verbose = TRUE
)

colMeans(MoEfitBayes$pi) # estimated pi
colMeans(dat$pi) # true pi

apply(MoEfitBayes$nu, 2, mean) # estimated nu
dat$nu # true nu


# run mixture model with EM algorithm
set.seed(123)
MoEfitEM <- moewishartX(
  S_list, X = matrix(rep(1, n), ncol = 1),
  method = "em",
  K = my_K_intial, #init_pi = init_pi, 
  niter = 10000, verbose = TRUE
)
# EM estimates
MoEfitEM$pi
colMeans(dat$pi) # true pi

MoEfitEM$nu
dat$nu # true nu

fit <- list(fitBayes = fitBayes, fitEM = fitEM, MoEfitBayes = MoEfitBayes, MoEfitEM = MoEfitEM)

save(fit, file = paste0(filename_save, seed.idx, ".RData"))


########################
## summarize results of mixture models
########################

n <- 200
filename_save <- paste0("/data/zhiz/moewishart/results/n", n, "_p2K3_sim")

n_sim <- 100
# for(i in 1:n_sim) {
#   load(paste0(filename_save, i, ".RData"))
#   cat("sim", i, "Bayesian pi=", colMeans(fit$fitBayes$pi), "; EM pi=", fit$fitEM$pi, "\n")
#   cat("sim", i, "Bayesian nu=", colMeans(fit$fitBayes$nu), "; EM nu=", fit$fitEM$nu, "\n")

#   cat(all(order(colMeans(fit$fitBayes$pi)) == order(colMeans(fit$fitBayes$nu))),
#     all(order(fit$fitEM$pi) == order(fit$fitEM$nu)), "\n")
# }


# library(moewishart)
# n <- 200 # subjects
# p <- 2
# K <- 3
# set.seed(1)
# dat <- simData(n, p,
#   pis = c(0.35, 0.40, 0.25),
#   nus = c(8, 12, 3),
#   Sigma = NULL # default Sigmas pre-defined in function 'simData()'
# )
# #S_list <- dat$S
# Sigma.true <- dat$Sigma_list
# pi.true <- dat$pi[1, ]
# nu.true <- dat$nu

# # pi errors
# pi.error <- matrix(nrow = n_sim, ncol = 6)
# nu.error <- matrix(nrow = n_sim, ncol = 6)
# Sigma.error <- matrix(nrow = n_sim, ncol = 6)
# colnames(pi.error) <- colnames(nu.error) <- colnames(Sigma.error) <- c("BayesMean1", "BayesMean2", "BMA1", "BMA2", "EM1", "EM2")

# for(i in 1:n_sim) {
#   load(paste0(filename_save, i, ".RData"))
#   # cat("sim", i, "Bayesian pi=", colMeans(fit$fitBayes$pi), "; EM pi=", fit$fitEM$pi, "\n")
#   # cat("sim", i, "Bayesian nu=", colMeans(fit$fitBayes$nu), "; EM nu=", fit$fitEM$nu, "\n")

#   TrueOrder <- order(nu.true)
#   BayesOrder <-  order(colMeans(fit$fitBayes$nu))
#   pi.error[i, 1] <- mean(abs(colMeans(fit$fitBayes$pi)[BayesOrder] - pi.true[TrueOrder]))
#   # pi.error[i, 2] <- sqrt(mean((colMeans(fit$fitBayes$pi)[BayesOrder] - pi.true[TrueOrder])^2))
#   pi.error[i, 2] <- mean((colMeans(fit$fitBayes$pi)[BayesOrder] - pi.true[TrueOrder])^2)
#   pi.error[i, 3] <- mean( apply(fit$fitBayes$pi[,BayesOrder], 1, function(xx) mean(abs(xx - pi.true[TrueOrder]))) )
#   # pi.error[i, 4] <- mean( apply(fit$fitBayes$pi[,BayesOrder], 1, function(xx) sqrt(mean((xx - pi.true[TrueOrder])^2))) )
#   pi.error[i, 4] <- mean( apply(fit$fitBayes$pi[,BayesOrder], 1, function(xx) mean((xx - pi.true[TrueOrder])^2)) )

#   nu.error[i, 1] <- mean(abs(colMeans(fit$fitBayes$nu)[BayesOrder] - nu.true[TrueOrder]))
#   # nu.error[i, 2] <- sqrt(mean((colMeans(fit$fitBayes$nu)[BayesOrder] - nu.true[TrueOrder])^2))
#   nu.error[i, 2] <- mean((colMeans(fit$fitBayes$nu)[BayesOrder] - nu.true[TrueOrder])^2)
#   nu.error[i, 3] <- mean( apply(fit$fitBayes$nu[,BayesOrder], 1, function(xx) mean(abs(xx - nu.true[TrueOrder]))) )
#   # nu.error[i, 4] <- mean( apply(fit$fitBayes$nu[,BayesOrder], 1, function(xx) sqrt(mean((xx - nu.true[TrueOrder])^2))) )
#   nu.error[i, 4] <- mean( apply(fit$fitBayes$nu[,BayesOrder], 1, function(xx) mean((xx - nu.true[TrueOrder])^2)) )

#   Sigma.bayes <- Reduce("+", fit$fitBayes$Sigma)/length(fit$fitBayes$Sigma)
#   Sigma.error[i, 1] <- mean( sapply(1:K, function(xx) sum(abs(Sigma.bayes[,,BayesOrder][,,xx] - Sigma.true[TrueOrder][[xx]])) ))
#   # Sigma.error[i, 2] <- sqrt(mean( sapply(1:K, function(xx) sum((Sigma.bayes[,,BayesOrder][,,xx] - Sigma.true[TrueOrder][[xx]])^2)) ))
#   # Sigma.error[i, 2] <- sqrt(mean( sapply(1:K, function(xx) (Sigma.bayes[,,BayesOrder][,,xx] - Sigma.true[TrueOrder][[xx]])^2) ))
#   Sigma.error[i, 2] <- mean( sapply(1:K, function(xx) {
#     s0 <- Sigma.true[TrueOrder][[xx]]
#     s <- Sigma.bayes[,,BayesOrder][,,xx]
#     sum((s[lower.tri(s, diag = TRUE)] - s0[lower.tri(s, diag = TRUE)])^2)
#   }) )
#   Sigma.error[i, 3] <- mean( unlist(lapply(fit$fitBayes$Sigma, function(Sigma_single){
#     mean( sapply(1:K, function(xx) sum(abs(Sigma_single[,,BayesOrder][,,xx] - Sigma.true[TrueOrder][[xx]])) ))
#   })) )
#   # Sigma.error[i, 4] <- mean( unlist(lapply(fit$fitBayes$Sigma, function(Sigma_single){
#   #   sqrt(mean( sapply(1:K, function(xx) (Sigma_single[,,BayesOrder][,,xx] - Sigma.true[TrueOrder][[xx]])^2) ))
#   # })) )
#   Sigma.error[i, 4] <- mean( unlist(lapply(fit$fitBayes$Sigma, function(Sigma_single){
#     sqrt(mean( sapply(1:K, function(xx) {
#     s0 <- Sigma.true[TrueOrder][[xx]]
#     s <- Sigma_single[,,BayesOrder][,,xx]
#     sum((s[lower.tri(s, diag = TRUE)] - s0[lower.tri(s, diag = TRUE)])^2)
#     } ) ))
#   })) )

#   EMorder <- order(fit$fitEM$nu)
#   pi.error[i, 5] <- mean(abs(fit$fitEM$pi[EMorder] - pi.true[TrueOrder]))
#   # pi.error[i, 6] <- sqrt(mean((fit$fitEM$pi[EMorder] - pi.true[TrueOrder])^2))
#   pi.error[i, 6] <- mean((fit$fitEM$pi[EMorder] - pi.true[TrueOrder])^2)
#   nu.error[i, 5] <- mean(abs(fit$fitEM$nu[EMorder] - nu.true[TrueOrder]))
#   # nu.error[i, 6] <- sqrt(mean((fit$fitEM$nu[EMorder] - nu.true[TrueOrder])^2))
#   nu.error[i, 6] <- mean((fit$fitEM$nu[EMorder] - nu.true[TrueOrder])^2)
#   Sigma.error[i, 5] <- mean( sapply(1:K, function(xx) sum(abs(fit$fitEM$Sigma[EMorder][[xx]] - Sigma.true[TrueOrder][[xx]]))) )
#   # Sigma.error[i, 6] <- sqrt(mean( sapply(1:K, function(xx) sum((fit$fitEM$Sigma[EMorder][[xx]] - Sigma.true[TrueOrder][[xx]])^2)) ))
#   Sigma.error[i, 6] <- mean( sapply(1:K, function(xx) sum((fit$fitEM$Sigma[EMorder][[xx]] - Sigma.true[TrueOrder][[xx]])^2)) )
# }

# # summary of means
# error.mean <- 
# rbind(
#   rbind(apply(pi.error, 2, mean)[c(1,3,5)], apply(pi.error, 2, mean)[c(2,4,6)]), 
#   rbind(apply(nu.error, 2, mean)[c(1,3,5)], apply(nu.error, 2, mean)[c(2,4,6)]),
#   rbind(apply(Sigma.error, 2, mean)[c(1,3,5)], apply(Sigma.error, 2, mean)[c(2,4,6)]) 
# )
# rownames(error.mean) <- c("pi-norm1", "pi-norm2", "nu-norm1", "nu-norm2", "Sigma-norm1", "Sigma-norm2")
# error.mean <- round(error.mean, digit = 3)

# # summary of sd
# error.sd <- 
# rbind(
#   rbind(apply(pi.error, 2, sd)[c(1,3,5)], apply(pi.error, 2, sd)[c(2,4,6)]), 
#   rbind(apply(nu.error, 2, sd)[c(1,3,5)], apply(nu.error, 2, sd)[c(2,4,6)]),
#   rbind(apply(Sigma.error, 2, sd)[c(1,3,5)], apply(Sigma.error, 2, sd)[c(2,4,6)]) 
# )
# error.sd <- round(error.sd, digit = 4)

# t( sapply(1:6, function(xx) as.vector(paste0(error.mean[xx, ], " (", error.sd[xx, ], ") "))) )



# ########################
# ## summarize results of MoE models (model misspecification)
# ########################

# # errors
# pi.error <- matrix(nrow = n_sim, ncol = 6)
# nu.error <- matrix(nrow = n_sim, ncol = 6)
# Sigma.error <- matrix(nrow = n_sim, ncol = 6)
# colnames(pi.error) <- colnames(nu.error) <- colnames(Sigma.error) <- c("BayesMean1", "BayesMean2", "BMA1", "BMA2", "EM1", "EM2")

# for(i in 1:n_sim) {
#   load(paste0(filename_save, i, ".RData"))

#   TrueOrder <- order(nu.true)
#   BayesOrder <-  order(colMeans(fit$MoEfitBayes$nu))
#   pi.error[i, 1] <- mean(abs(colMeans(fit$MoEfitBayes$pi)[BayesOrder] - pi.true[TrueOrder]))
#   # pi.error[i, 2] <- sqrt(mean((colMeans(fit$MoEfitBayes$pi)[BayesOrder] - pi.true[TrueOrder])^2))
#   pi.error[i, 2] <- mean((colMeans(fit$MoEfitBayes$pi)[BayesOrder] - pi.true[TrueOrder])^2)
#   pi.error[i, 3] <- mean( apply(fit$MoEfitBayes$pi[,BayesOrder], 1, function(xx) mean(abs(xx - pi.true[TrueOrder]))) )
#   # pi.error[i, 4] <- mean( apply(fit$MoEfitBayes$pi[,BayesOrder], 1, function(xx) sqrt(mean((xx - pi.true[TrueOrder])^2))) )
#   pi.error[i, 4] <- mean( apply(fit$MoEfitBayes$pi[,BayesOrder], 1, function(xx) mean((xx - pi.true[TrueOrder])^2)) )

#   nu.error[i, 1] <- mean(abs(colMeans(fit$MoEfitBayes$nu)[BayesOrder] - nu.true[TrueOrder]))
#   # nu.error[i, 2] <- sqrt(mean((colMeans(fit$MoEfitBayes$nu)[BayesOrder] - nu.true[TrueOrder])^2))
#   nu.error[i, 2] <- mean((colMeans(fit$MoEfitBayes$nu)[BayesOrder] - nu.true[TrueOrder])^2)
#   nu.error[i, 3] <- mean( apply(fit$MoEfitBayes$nu[,BayesOrder], 1, function(xx) mean(abs(xx - nu.true[TrueOrder]))) )
#   # nu.error[i, 4] <- mean( apply(fit$MoEfitBayes$nu[,BayesOrder], 1, function(xx) sqrt(mean((xx - nu.true[TrueOrder])^2))) )
#   nu.error[i, 4] <- mean( apply(fit$MoEfitBayes$nu[,BayesOrder], 1, function(xx) mean((xx - nu.true[TrueOrder])^2)) )

#   Sigma.bayes <- Reduce("+", fit$MoEfitBayes$Sigma)/length(fit$MoEfitBayes$Sigma)
#   Sigma.error[i, 1] <- mean( sapply(1:K, function(xx) sum(abs(Sigma.bayes[,,BayesOrder][,,xx] - Sigma.true[TrueOrder][[xx]])) ))
#   # Sigma.error[i, 2] <- sqrt(mean( sapply(1:K, function(xx) sum((Sigma.bayes[,,BayesOrder][,,xx] - Sigma.true[TrueOrder][[xx]])^2)) ))
#   # Sigma.error[i, 2] <- sqrt(mean( sapply(1:K, function(xx) (Sigma.bayes[,,BayesOrder][,,xx] - Sigma.true[TrueOrder][[xx]])^2) ))
#   Sigma.error[i, 2] <- mean( sapply(1:K, function(xx) {
#     s0 <- Sigma.true[TrueOrder][[xx]]
#     s <- Sigma.bayes[,,BayesOrder][,,xx]
#     sum((s[lower.tri(s, diag = TRUE)] - s0[lower.tri(s, diag = TRUE)])^2)
#   }) )
#   Sigma.error[i, 3] <- mean( unlist(lapply(fit$MoEfitBayes$Sigma, function(Sigma_single){
#     mean( sapply(1:K, function(xx) sum(abs(Sigma_single[,,BayesOrder][,,xx] - Sigma.true[TrueOrder][[xx]])) ))
#   })) )
#   # Sigma.error[i, 4] <- mean( unlist(lapply(fit$MoEfitBayes$Sigma, function(Sigma_single){
#   #   sqrt(mean( sapply(1:K, function(xx) (Sigma_single[,,BayesOrder][,,xx] - Sigma.true[TrueOrder][[xx]])^2) ))
#   # })) )
#   Sigma.error[i, 4] <- mean( unlist(lapply(fit$MoEfitBayes$Sigma, function(Sigma_single){
#     sqrt(mean( sapply(1:K, function(xx) {
#     s0 <- Sigma.true[TrueOrder][[xx]]
#     s <- Sigma_single[,,BayesOrder][,,xx]
#     sum((s[lower.tri(s, diag = TRUE)] - s0[lower.tri(s, diag = TRUE)])^2)
#     } ) ))
#   })) )

#   EMorder <- order(fit$MoEfitEM$nu)
#   pi.error[i, 5] <- mean(abs(fit$MoEfitEM$gamma[EMorder] - pi.true[TrueOrder]))
#   # pi.error[i, 6] <- sqrt(mean((fit$MoEfitEM$pi[EMorder] - pi.true[TrueOrder])^2))
#   pi.error[i, 6] <- mean((fit$MoEfitEM$gamma[EMorder] - pi.true[TrueOrder])^2)
#   nu.error[i, 5] <- mean(abs(fit$MoEfitEM$nu[EMorder] - nu.true[TrueOrder]))
#   # nu.error[i, 6] <- sqrt(mean((fit$MoEfitEM$nu[EMorder] - nu.true[TrueOrder])^2))
#   nu.error[i, 6] <- mean((fit$MoEfitEM$nu[EMorder] - nu.true[TrueOrder])^2)
#   Sigma.error[i, 5] <- mean( sapply(1:K, function(xx) sum(abs(fit$MoEfitEM$Sigma[EMorder][[xx]] - Sigma.true[TrueOrder][[xx]]))) )
#   # Sigma.error[i, 6] <- sqrt(mean( sapply(1:K, function(xx) sum((fit$MoEfitEM$Sigma[EMorder][[xx]] - Sigma.true[TrueOrder][[xx]])^2)) ))
#   Sigma.error[i, 6] <- mean( sapply(1:K, function(xx) sum((fit$MoEfitEM$Sigma[EMorder][[xx]] - Sigma.true[TrueOrder][[xx]])^2)) )
# }

# # summary of means
# error.mean <- 
# rbind(
#   rbind(apply(pi.error, 2, mean)[c(1,3,5)], apply(pi.error, 2, mean)[c(2,4,6)]), 
#   rbind(apply(nu.error, 2, mean)[c(1,3,5)], apply(nu.error, 2, mean)[c(2,4,6)]),
#   rbind(apply(Sigma.error, 2, mean)[c(1,3,5)], apply(Sigma.error, 2, mean)[c(2,4,6)]) 
# )
# rownames(error.mean) <- c("pi-norm1", "pi-norm2", "nu-norm1", "nu-norm2", "Sigma-norm1", "Sigma-norm2")
# error.mean <- round(error.mean, digit = 3)

# # summary of sd
# error.sd <- 
# rbind(
#   rbind(apply(pi.error, 2, sd)[c(1,3,5)], apply(pi.error, 2, sd)[c(2,4,6)]), 
#   rbind(apply(nu.error, 2, sd)[c(1,3,5)], apply(nu.error, 2, sd)[c(2,4,6)]),
#   rbind(apply(Sigma.error, 2, sd)[c(1,3,5)], apply(Sigma.error, 2, sd)[c(2,4,6)]) 
# )
# error.sd <- round(error.sd, digit = 4)

# t( sapply(1:6, function(xx) as.vector(paste0(error.mean[xx, ], " (", error.sd[xx, ], ") "))) )


