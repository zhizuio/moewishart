rm(list = ls())

n <- 1000 # subjects
p <- 2

filename_save <- paste0("/data/zhiz/moewishart/results/n", n, "_p2K3_sim")
args <- commandArgs(trailingOnly = TRUE)
seed.idx <- as.numeric(args)

con <- file(paste0(filename_save, seed.idx,".log"))
sink(con, append=TRUE)
sink(con, append=TRUE, type="message")


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

# run Bayesian moewishart model
set.seed(123)
fitBayes <- moewishart(S_list,
  K = my_K_intial, init_pi = init_pi, 
  mh_sigma = 0.2,
  niter = 10000, burnin = 1000, thin = 1, verbose = TRUE
)

colMeans(fitBayes$pi) # estimated pi
colMeans(dat$pi) # true pi

apply(fitBayes$nu, 2, mean) # estimated nu
dat$nu # true nu


# run the model with EM algorithm
set.seed(123)
fitEM <- moewishart(S_list,
  method = "em",
  K = my_K_intial, init_pi = init_pi, 
  niter = 10000, verbose = TRUE
)
# EM estimates
fitEM$pi
colMeans(dat$pi) # true pi

fitEM$nu
dat$nu # true nu

fit <- list(fitBayes = fitBayes, fitEM = fitEM)

save(fit, file = paste0(filename_save, seed.idx, ".RData"))

########################
## summarize results
########################

n <- 1000
filename_save <- paste0("/data/zhiz/moewishart/results/n", n, "_p2K3_sim")

n_sim <- 100
## for(i in 1:n_sim) {
##   load(paste0(filename_save, i, ".RData"))
##   cat("sim", i, "Bayesian pi=", colMeans(fit$fitBayes$pi), "; EM pi=", fit$fitEM$pi, "\n")
##   cat("sim", i, "Bayesian nu=", colMeans(fit$fitBayes$nu), "; EM nu=", fit$fitEM$nu, "\n")

##   cat(all(order(colMeans(fit$fitBayes$pi)) == order(colMeans(fit$fitBayes$nu))),
##     all(order(fit$fitEM$pi) == order(fit$fitEM$nu)), "\n")
## }


# library(moewishart)
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
# pi.error <- matrix(nrow = n_sim, ncol = 3)
# nu.error <- matrix(nrow = n_sim, ncol = 3)
# Sigma.error <- matrix(nrow = n_sim, ncol = 3)
# colnames(pi.error) <- colnames(nu.error) <- colnames(Sigma.error) <- c("BayesMean", "BMA", "EM")

# for(i in 1:n_sim) {
#   load(paste0(filename_save, i, ".RData"))
#   # cat("sim", i, "Bayesian pi=", colMeans(fit$fitBayes$pi), "; EM pi=", fit$fitEM$pi, "\n")
#   # cat("sim", i, "Bayesian nu=", colMeans(fit$fitBayes$nu), "; EM nu=", fit$fitEM$nu, "\n")

#   TrueOrder <- order(nu.true)
#   BayesOrder <-  order(colMeans(fit$fitBayes$nu))
#   pi.error[i, 1] <- sqrt(mean((colMeans(fit$fitBayes$pi)[BayesOrder] - pi.true[TrueOrder])^2))
#   pi.error[i, 2] <- mean( apply(fit$fitBayes$pi[,BayesOrder], 1, function(xx) sqrt(mean((xx - pi.true[TrueOrder])^2))) )
#   nu.error[i, 1] <- sqrt(mean((colMeans(fit$fitBayes$nu)[BayesOrder] - nu.true[TrueOrder])^2))
#   nu.error[i, 2] <- mean( apply(fit$fitBayes$nu[,BayesOrder], 1, function(xx) sqrt(mean((xx - nu.true[TrueOrder])^2))) )

#   Sigma.bayes <- Reduce("+", fit$fitBayes$Sigma)/length(fit$fitBayes$Sigma)
#   Sigma.error[i, 1] <- sqrt(mean( sapply(1:K, function(xx) sum((Sigma.bayes[,,BayesOrder][,,xx] - Sigma.true[TrueOrder][[xx]])^2)) ))
#   Sigma.error[i, 2] <- mean( unlist(lapply(fit$fitBayes$Sigma, function(Sigma_single){
#     sqrt(mean( sapply(1:K, function(xx) sum((Sigma_single[,,BayesOrder][,,xx] - Sigma.true[TrueOrder][[xx]])^2)) ))
#   })) )

#   EMorder <- order(fit$fitEM$nu)
#   pi.error[i, 3] <- sqrt(mean((fit$fitEM$pi[EMorder] - pi.true[TrueOrder])^2))
#   nu.error[i, 3] <- sqrt(mean((fit$fitEM$nu[EMorder] - nu.true[TrueOrder])^2))
#   Sigma.error[i, 3] <- sqrt(mean( sapply(1:K, function(xx) sum((fit$fitEM$Sigma[EMorder][[xx]] - Sigma.true[TrueOrder][[xx]])^2)) ))
# }

# apply(pi.error, 2, mean)
# apply(pi.error, 2, sd)

# apply(nu.error, 2, mean)
# apply(nu.error, 2, sd)

# apply(Sigma.error, 2, mean)
# apply(Sigma.error, 2, sd)
