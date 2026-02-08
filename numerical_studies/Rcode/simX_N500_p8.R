##========================================================================================
## This file contains code for 100 simulations
## Simulated data: Generated from "MoE Model & n=500 & p=8"
## Working models: Bayes, EM, Bayes-MoE, EM-MoE
## 
## Output 1: Each simulation fitted by the 4 methods -> one saved file 
## Output 2: Compute all errors of estimates from the 100 simulations and save them in one file
##========================================================================================

rm(list = ls())

# user-defined directory to save simulation results as data files
file_directory <- "../moewishart/"
filename_save <- paste0(file_directory, "n500_p8_sim")

# simulation settings
n <- 500 # subjects
p <- 8 # covariance dimension
K <- 3 # number of mixture components
n_sim <- 100 # number of simulations


make_pd_AR <- function(p, rho = 0.2) {
  Sigma <- outer(1:p, 1:p, function(i, j)rho^abs(i - j))
}
set.seed(123) # fix coefficients of underlying MoE model
Xq <- 3; K = 3
betas <- matrix(runif(Xq * K, -2, 2), nrow = Xq, ncol = K)
betas[, K] <- 0

library(moewishart)

for (seed.idx in 1:n_sim) {
  
  # simulate data from the underlying true model: "MoE model"
  set.seed(seed.idx)
  dat <- simData(n, p,
                 Xq = 3, K =3, betas = betas,
                 pis = c(0.35, 0.40, 0.25),
                 nus = c(9, 20, 14),
                 Sigma = list(  make_pd_AR(p,rho=0.5), make_pd_AR(p,rho=0.2), make_pd_AR(p,rho=.8) )
  )
  S_list <- dat$S
  Sigma_list <- dat$Sigma_list
  my_K_intial <- 3
    
  # run Bayesian mixture model
  set.seed(123)
  fitBayes <- moewishart(S_list,
    K = my_K_intial,  
    nu_prior_a = 2, nu_prior_b = 0.5,
    mh_sigma = 0.04,
    niter = 20000, burnin = 5000, thin = 1, verbose = TRUE
  )

  # run mixture model with EM algorithm
  set.seed(123)
  fitEM <- moewishart(S_list,
    method = "em", #cpp = TRUE,
    K = my_K_intial,  
    niter = 10000, verbose = TRUE,
    n_restarts = 20, restart_iters = 40, tol = 1e-8
  )

  # run Bayesian MoE model
  set.seed(123)
  MoEfitBayes <- moewishartX(
    S_list, X = cbind(1, dat$X),
    K = my_K_intial,  
    nu_prior_a = 2, nu_prior_b = 0.5,
    mh_sigma = 0.04, mh_beta = 0.13,
    niter = 20000, burnin = 5000, thin = 1, verbose = TRUE
  )

  # run MoE model with EM algorithm
  set.seed(123)
  MoEfitEM <- moewishartX(
    S_list, X = cbind(1, dat$X),
    method = "em",
    K = my_K_intial,  
    niter = 10000, verbose = TRUE
  )
  
  # put together the results of all 4 methods
  fit <- list(fitBayes = fitBayes, fitEM = fitEM, MoEfitBayes = MoEfitBayes, MoEfitEM = MoEfitEM)
  
  save(fit, file = paste0(filename_save, seed.idx, ".RData"))
}

########################
## summarize results
########################

burnin <- 5000

# true parameters
pi.true <- dat$pi[1, ]
nu.true <- dat$nu
Sigma.true <- dat$Sigma_list
beta.true <- dat$betas

# compute error of each estimate by methods 'Bayes' and 'EM'
# please be sure that the following R code file "moewishartX_errors.R" is ready
source("moewishart_errors.R")

# compute error of each estimate by methods 'Bayes-MoE' and 'EM-MoE'
# please be sure that the following R code file "moewishartX_errors.R" is ready
source("moewishartX_errors.R")

# save results as a data frame
datSimN500 <- data.frame(n = NULL, method = NULL, estimator = NULL, error = NULL )
datSimN500 <- rbind(datSimN500, data.frame(n = "n=500", method = "Bayes0", estimator = "pi-norm1", error = pi.error[,1] )) # posterior mean
datSimN500 <- rbind(datSimN500, data.frame(n = "n=500", method = "EM", estimator = "pi-norm1", error = pi.error[,3] ))
datSimN500 <- rbind(datSimN500, data.frame(n = "n=500", method = "Bayes0.X", estimator = "pi-norm1", error = pi.error2[,1] ))
datSimN500 <- rbind(datSimN500, data.frame(n = "n=500", method = "EM.X", estimator = "pi-norm1", error = pi.error2[,3] ))
datSimN500 <- rbind(datSimN500, data.frame(n = "n=500", method = "Bayes0", estimator = "nu-norm1", error = nu.error[,1] ))
datSimN500 <- rbind(datSimN500, data.frame(n = "n=500", method = "EM", estimator = "nu-norm1", error = nu.error[,3] ))
datSimN500 <- rbind(datSimN500, data.frame(n = "n=500", method = "Bayes0.X", estimator = "nu-norm1", error = nu.error2[,1] ))
datSimN500 <- rbind(datSimN500, data.frame(n = "n=500", method = "EM.X", estimator = "nu-norm1", error = nu.error2[,3] ))
datSimN500 <- rbind(datSimN500, data.frame(n = "n=500", method = "Bayes0", estimator = "Sigma-norm1", error = Sigma.error[,1] ))
datSimN500 <- rbind(datSimN500, data.frame(n = "n=500", method = "EM", estimator = "Sigma-norm1", error = Sigma.error[,3] ))
datSimN500 <- rbind(datSimN500, data.frame(n = "n=500", method = "Bayes0.X", estimator = "Sigma-norm1", error = Sigma.error2[,1] ))
datSimN500 <- rbind(datSimN500, data.frame(n = "n=500", method = "EM.X", estimator = "Sigma-norm1", error = Sigma.error2[,3] ))
datSimN500 <- rbind(datSimN500, data.frame(n = "n=500", method = "Bayes0", estimator = "Sigma-normF", error = Sigma.error[,2] ))
datSimN500 <- rbind(datSimN500, data.frame(n = "n=500", method = "EM", estimator = "Sigma-normF", error = Sigma.error[,4] ))
datSimN500 <- rbind(datSimN500, data.frame(n = "n=500", method = "Bayes0.X", estimator = "Sigma-normF", error = Sigma.error2[,2] ))
datSimN500 <- rbind(datSimN500, data.frame(n = "n=500", method = "EM.X", estimator = "Sigma-normF", error = Sigma.error2[,4] ))

datSimN500 <- rbind(datSimN500, data.frame(n = "n=500", method = "Bayes0.X", estimator = "Beta-norm1", error = beta.error2[,1] ))
datSimN500 <- rbind(datSimN500, data.frame(n = "n=500", method = "Bayes0.X", estimator = "Beta-normF", error = beta.error2[,2] ))
datSimN500 <- rbind(datSimN500, data.frame(n = "n=500", method = "EM.X", estimator = "Beta-norm1", error = beta.error2[,3] ))
datSimN500 <- rbind(datSimN500, data.frame(n = "n=500", method = "EM.X", estimator = "Beta-normF", error = beta.error2[,4] ))

save(datSimN500, file = paste0(file_directory, "/datSimX_N500.RData"))
