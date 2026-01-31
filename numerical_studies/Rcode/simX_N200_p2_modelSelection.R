##========================================================================================
## This file contains code for 100 simulations for model selection
## Simulated data: Generated from "MoE Model & n=200 & p=2"
## Working models: Bayes, EM, Bayes-MoE, EM-MoE with K={2,3,4,5,6}
## 
## Output: AIC, BIC, ICL, elpd_loo
##========================================================================================

rm(list = ls())

# user-defined directory to save simulation results as data files
file_directory <- "../moewishartX/"
setwd(file_directory)

# summarize information criteria
datIC <- data.frame(model = NULL, K = NULL, simIdx = NULL, elpd = NULL, icl = NULL, icl2 = NULL, aic = NULL, bic = NULL)
niter <- 20000; burnin <- 5000; thin <- 1
digits0 <- 0

# simulation settings
n <- 200 # subjects
p <- 2 # covariance dimension
q <- 0 # covariates 
K <- 3 # number of mixture components
set.seed(123) # fix coefficients of underlying MoE model
q <- 3; K = 3
betas <- matrix(runif(q * K, -2, 2), nrow = q, ncol = K)
betas[, K] <- 0

n_sim <- 100 # number of simulations

library(moewishart)

for (my_K_intial in 2:6) {

  filename_save <- paste0(file_directory, "n", n, "_p", p, "K", my_K_intial, "_sim")

  for (seed.idx in 1:n_sim) {
    
    # simulate data from the underlying true model: "mixture of covariance matrix model"
    set.seed(seed.idx)
    dat <- simData(n, p,
      Xq = 3, K = 3, betas = betas,
      pis = c(0.35, 0.40, 0.25),
      nus = c(8, 12, 3),
      Sigma = NULL # default Sigmas pre-defined in function 'simData()'
    )
    S_list <- dat$S
    Sigma_list <- dat$Sigma_list
    
    # run Bayesian mixture model
    set.seed(123)
    fitBayes <- moewishart(S_list,
      K = my_K_intial,  
      nu_prior_a = 2, nu_prior_b = 0.5,
      mh_sigma = 0.2,
      niter = 20000, burnin = 5000, thin = 1, verbose = TRUE
    )
    a <- loo::loo(fitBayes$loglik_individual[((1 + burnin):niter) / thin, ])
    fitBayes$elpd <- a$estimates["elpd_loo",]

    # run mixture model with EM algorithm
    set.seed(123)
    fitEM <- moewishart(S_list,
      method = "em", #cpp = TRUE,
      K = my_K_intial,  
      niter = 10000, verbose = TRUE,
      n_restarts = 5, restart_iters = 40, tol = 1e-8
    )

    # run Bayesian MoE model
    set.seed(123)
    MoEfitBayes <- moewishartX(
      S_list, X = cbind(1, dat$X),
      K = my_K_intial,  
      nu_prior_a = 2, nu_prior_b = 0.5,
      mh_sigma = 0.2, mh_beta = 0.25,
      niter = 20000, burnin = 5000, thin = 1, verbose = TRUE
    )
    a <- loo::loo(MoEfitBayes$loglik_individual[((1 + burnin):niter) / thin, ])
    MoEfitBayes$elpd <- a$estimates["elpd_loo",]

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

    # load(paste0(filename_save, seed.idx, ".RData"))
    

    ########################
    ## compute information criteria
    ########################

    # Bayes MM
    elpd <- fit$fitBayes$elpd[1]

    bic <- -2 * sum(fit$fitBayes$loglik[-c(1:burnin)]) + log(n) * (my_K_intial + p*(p+1)/2 + (my_K_intial-1)*(q+1)) 
    r_ik <- fit$fitBayes$pi_ik[-c(1:burnin), , ]
    r_ik <- apply(r_ik, c(2, 3), mean)
    icl <- bic - 2 * sum(r_ik * log(r_ik), na.rm = TRUE)
    # icl <- round(icl, digits = digits0)
    c_ik <- t( apply(r_ik, 1, function(xx) as.numeric(xx==max(xx))) )
    icl2 <- bic - 2 * sum(c_ik * log(r_ik), na.rm = TRUE)
    
    datIC <- rbind(datIC, data.frame(model = "Bayes", K = my_K_intial, simIdx = seed.idx, elpd = elpd, icl = icl, icl2 = icl2, aic = NA, bic = NA))

    # Bayes MoE
    elpd <- fit$MoEfitBayes$elpd[1]

    bic <- -2 * sum(fit$MoEfitBayes$loglik[-c(1:burnin)]) + log(n) * (my_K_intial + p*(p+1)/2 + (my_K_intial-1)*(q+1)) 
    r_ik <- fit$MoEfitBayes$pi_ik[-c(1:burnin), , ]
    r_ik <- apply(r_ik, c(2, 3), mean)
    icl <- bic - 2 * sum(r_ik * log(r_ik), na.rm = TRUE)
    # icl <- round(icl, digits = digits0)
    c_ik <- t( apply(r_ik, 1, function(xx) as.numeric(xx==max(xx))) )
    icl2 <- bic - 2 * sum(c_ik * log(r_ik), na.rm = TRUE)
    
    datIC <- rbind(datIC, data.frame(model = "Bayes-MoE", K = my_K_intial, simIdx = seed.idx, elpd = elpd, icl = icl, icl2 = icl2, aic = NA, bic = NA))


    # EM MM
    aic <- -2 * sum(fit$fitEM$loglik) + 2 * (my_K_intial * (p*(p+1)/2+1) + (my_K_intial-1) * 1)
    bic <- -2 * sum(fit$fitEM$loglik) + log(n) * (my_K_intial * (p*(p+1)/2+1) + (my_K_intial-1) * q)
    r_ik <- fit$fitEM$tau
    icl <- bic - 2 * sum(r_ik * log(r_ik), na.rm = TRUE)
    # aic <- round(aic, digits = digits0)
    # bic <- round(bic, digits = digits0)
    # icl <- round(icl, digits = digits0)
    c_ik <- t( apply(r_ik, 1, function(xx) as.numeric(xx==max(xx))) )
    icl2 <- bic - 2 * sum(c_ik * log(r_ik), na.rm = TRUE)
    
    datIC <- rbind(datIC, data.frame(model = "EM", K = my_K_intial, simIdx = seed.idx, elpd = NA, icl = icl, icl2 = icl2, aic = aic, bic = bic))

    # EM MoE
    aic <- -2 * sum(fit$MoEfitEM$loglik) + 2 * (my_K_intial * (p*(p+1)/2+1) + (my_K_intial-1) * 1)
    bic <- -2 * sum(fit$MoEfitEM$loglik) + log(n) * (my_K_intial * (p*(p+1)/2+1) + (my_K_intial-1) * q)
    r_ik <- fit$MoEfitEM$gamma
    icl <- bic - 2 * sum(r_ik * log(r_ik), na.rm = TRUE)
    # aic <- round(aic, digits = digits0)
    # bic <- round(bic, digits = digits0)
    # icl <- round(icl, digits = digits0)
    c_ik <- t( apply(r_ik, 1, function(xx) as.numeric(xx==max(xx))) )
    icl2 <- bic - 2 * sum(c_ik * log(r_ik), na.rm = TRUE)
    
    datIC <- rbind(datIC, data.frame(model = "EM-MoE", K = my_K_intial, simIdx = seed.idx, elpd = NA, icl = icl, icl2 = icl2, aic = aic, bic = bic))

  }
}


IC_mean <- data.frame("model" = NULL, K = NULL, elpd = NULL, icl = NULL, icl2 = NULL, aic = NULL, bic = NULL)
IC_sd <- data.frame("model" = NULL, K = NULL, elpd = NULL, icl = NULL, icl2 = NULL, aic = NULL, bic = NULL)
for (model in unique(datIC$model)){
  for (k in 2:6) {
    tmp <- datIC[datIC$model == model & datIC$K == k, 4:8]
    IC_mean <- rbind(IC_mean, data.frame("model" = model, K = k, t(data.frame(apply(tmp, 2, mean)))))
    IC_sd <- rbind(IC_sd, data.frame("model" = model, K = k, t(data.frame(apply(tmp, 2, sd)))))
  }
}
rownames(IC_mean) <- NULL
rownames(IC_mean) <- NULL

IC_mean_sd <- IC_mean
IC_mean_sd[, 3:7] <- matrix(paste0(as.matrix(round(IC_mean[, 3:7],2)), " (", as.matrix(round(IC_sd[, 3:7],2)), ")"), nrow=nrow(IC_mean))

datIC_all <- list(datIC = datIC, IC_mean = IC_mean, IC_sd = IC_sd)
save(datIC_all, file = paste0(file_directory, "n200_p", p, "datIC.RData"))

