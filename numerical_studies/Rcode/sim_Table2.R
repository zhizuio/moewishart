##========================================================================================
## This file contains code to produce Table 2
## 
## Input: Data files from R code 
##        'simX_N200_p2.R', 'simX_N500_p2.R', 'simX_N1000_p2.R'
##        'simX_N200_p8.R', 'simX_N500_p8.R', 'simX_N1000_p8.R'
## Output: Table 1
##========================================================================================

rm(list = ls())

# specify folder directory including the results from R code 'simX_N200_p2.R', 'simX_N500_p2.R' etc.
file_directory <- "../moewishart/"

##=================================
### True model: Mixture Model & p=2
##=================================

p <- 2
burnin <- 1000
ESS.MoE <- ESS.MM <- data.frame(p = NULL, n = NULL, ESS = NULL)
ESS.MoE.sd <- ESS.MM.sd <- data.frame(p = NULL, n = NULL, ESS = NULL)

n_all <- c(200, 500, 1000) # subjects
for (n in n_all) {
  
  filename_save <- paste0(file_directory, "n", n, "_p", p, "_sim")
  
  n_sim <- 100
  
  library(coda)
  library(moewishart)
  K <- 3
  make_pd_AR <- function(p, rho = 0.2) {
    Sigma <- outer(1:p, 1:p, function(i, j)rho^abs(i - j))
  }
  set.seed(1)
  dat <- simData(n, p,
                 pis = c(0.35, 0.40, 0.25),
                 nus = c(9, 20, 14),
                 Sigma = list(  make_pd_AR(p,rho=0.5), make_pd_AR(p,rho=0.2), make_pd_AR(p,rho=.8) )
  )
  Sigma.true <- dat$Sigma_list
  pi.true <- dat$pi[1, ]
  nu.true <- dat$nu
  
  ess.MM.pi <- ess.MoE.pi <- matrix(nrow = n_sim, ncol = 3)
  ess.MM.nu <- ess.MoE.nu <- matrix(nrow = n_sim, ncol = 3)
  ess.MM.Sigma <- ess.MoE.Sigma <- matrix(nrow = n_sim, ncol = 3) # only show 3 entries, i.e. Sigma11_1, Sigma11_2, Sigma11_3
  ess.MM.beta <- ess.MoE.beta <- matrix(nrow = n_sim, ncol = 2) # only show 2 entries, i.e. beta_11, beta_12
  
  for(i in 1:n_sim) {
    load(paste0(filename_save, i, ".RData"))
    
    ## ESS for MCMC from mixture model
    BayesOrder <-  order(colMeans(fit$fitBayes$nu[-c(1:burnin),]))
    
    ess.MM.pi[i, ] <- coda::effectiveSize(as.mcmc(fit$fitBayes$pi[-c(1:burnin),BayesOrder]))
    ess.MM.nu[i, ] <- coda::effectiveSize(as.mcmc(fit$fitBayes$nu[-c(1:burnin),BayesOrder]))
    
    ## extract the first Sigma's entry from each of the 3 components
    Sigmas <- cbind(unlist(lapply(fit$fitBayes$Sigma[-c(1:burnin)], function(xx) xx[1, 1, BayesOrder[1]])), 
                    unlist(lapply(fit$fitBayes$Sigma[-c(1:burnin)], function(xx) xx[1, 1, BayesOrder[2]])), 
                    unlist(lapply(fit$fitBayes$Sigma[-c(1:burnin)], function(xx) xx[1, 1, BayesOrder[3]])) )
    ess.MM.Sigma[i, ] <- coda::effectiveSize(as.mcmc(Sigmas))
    
    
    ## ESS for MCMC from MoE model
    BayesOrder <-  order(colMeans(fit$MoEfitBayes$nu[-c(1:burnin),]))
    
    ess.MoE.nu[i, ] <- coda::effectiveSize(as.mcmc(fit$MoEfitBayes$nu[-c(1:burnin),BayesOrder]))
    
    ## extract the first Sigma's entry from each of the 3 components
    Sigmas <- cbind(unlist(lapply(fit$MoEfitBayes$Sigma[-c(1:burnin)], function(xx) xx[1, 1, BayesOrder[1]])), 
                    unlist(lapply(fit$MoEfitBayes$Sigma[-c(1:burnin)], function(xx) xx[1, 1, BayesOrder[2]])), 
                    unlist(lapply(fit$MoEfitBayes$Sigma[-c(1:burnin)], function(xx) xx[1, 1, BayesOrder[3]])) )
    ess.MoE.Sigma[i, ] <- coda::effectiveSize(as.mcmc(Sigmas))
    
    if( dim(fit$MoEfitBayes$Beta_samples)[2] > 1 ) {
      # extract the first coefficients from two components; always the last dimension's value 0
      Betas <- fit$MoEfitBayes$Beta_samples[-c(1:burnin), 1, 1:2] 
      ess.MoE.beta[i, ] <- coda::effectiveSize(as.mcmc(Betas))
    }
    
  }
  
  ### calculate Mean and Standard Deviations
  value <- paste(c(
    paste(paste0(round(colMeans(ess.MM.nu), 0), " (", round(apply(ess.MM.nu, 2, sd), 0), ")")[1:2], collapse = " & "),
    paste(paste0(round(colMeans(ess.MM.Sigma), 0), " (", round(apply(ess.MM.Sigma, 2, sd), 0), ")")[1:2], collapse = " & "),
    paste(paste0(round(colMeans(ess.MM.pi), 0), " (", round(apply(ess.MM.pi, 2, sd), 0), ")")[1:2], collapse = " & ")
  ), collapse = " && ")
  ESS.MM.sd <- rbind(ESS.MM.sd, data.frame(p = p, n = paste0("n=",n), ESS = value))
  
  value <- paste(c(
    paste(paste0(round(colMeans(ess.MoE.nu), 0), " (", round(apply(ess.MoE.nu, 2, sd), 0), ")")[1:2], collapse = " & "),
    paste(paste0(round(colMeans(ess.MoE.Sigma), 0), " (", round(apply(ess.MoE.Sigma, 2, sd), 0), ")")[1:2], collapse = " & ")
  ), collapse = " && ")
  ESS.MoE.sd <- rbind(ESS.MoE.sd, data.frame(p = p, n = paste0("n=",n), ESS = value))
  
  
}


##=================================
### True model: Mixture Model & p=8
##=================================

p <- 8
burnin <- 5000
#ESS.MoE <- ESS.MM <- data.frame(p = NULL, n = NULL, ESS = NULL)

n_all <- c(200, 500, 1000) # subjects
for (n in n_all) {
  
  filename_save <- paste0(file_directory, "n", n, "_p", p, "_sim")
  
  n_sim <- 100
  
  library(coda)
  library(moewishart)
  p <- 8
  K <- 3
  make_pd_AR <- function(p, rho = 0.2) {
    Sigma <- outer(1:p, 1:p, function(i, j)rho^abs(i - j))
  }
  set.seed(1)
  dat <- simData(n, p,
                 pis = c(0.35, 0.40, 0.25),
                 nus = c(9, 20, 14),
                 Sigma = list(  make_pd_AR(p,rho=0.5), make_pd_AR(p,rho=0.2), make_pd_AR(p,rho=.8) )
  )
  Sigma.true <- dat$Sigma_list
  pi.true <- dat$pi[1, ]
  nu.true <- dat$nu
  
  ess.MM.pi <- ess.MoE.pi <- matrix(nrow = n_sim, ncol = 3)
  ess.MM.nu <- ess.MoE.nu <- matrix(nrow = n_sim, ncol = 3)
  ess.MM.Sigma <- ess.MoE.Sigma <- matrix(nrow = n_sim, ncol = 3) # only show 3 entries, i.e. Sigma11_1, Sigma11_2, Sigma11_3
  ess.MM.beta <- ess.MoE.beta <- matrix(nrow = n_sim, ncol = 2) # only show 2 entries, i.e. beta_11, beta_12
  
  for(i in 1:n_sim) {
    load(paste0(filename_save, i, ".RData"))
    
    ## ESS for MCMC from mixture model
    BayesOrder <-  order(colMeans(fit$fitBayes$nu[-c(1:burnin),]))
    
    ess.MM.pi[i, ] <- coda::effectiveSize(as.mcmc(fit$fitBayes$pi[-c(1:burnin),BayesOrder]))
    ess.MM.nu[i, ] <- coda::effectiveSize(as.mcmc(fit$fitBayes$nu[-c(1:burnin),BayesOrder]))
    
    ## extract the first Sigma's entry from each of the 3 components
    Sigmas <- cbind(unlist(lapply(fit$fitBayes$Sigma[-c(1:burnin)], function(xx) xx[1, 1, BayesOrder[1]])), 
                    unlist(lapply(fit$fitBayes$Sigma[-c(1:burnin)], function(xx) xx[1, 1, BayesOrder[2]])), 
                    unlist(lapply(fit$fitBayes$Sigma[-c(1:burnin)], function(xx) xx[1, 1, BayesOrder[3]])) )
    ess.MM.Sigma[i, ] <- coda::effectiveSize(as.mcmc(Sigmas))
    
    
    ## ESS for MCMC from MoE model
    BayesOrder <-  order(colMeans(fit$MoEfitBayes$nu[-c(1:burnin),]))
    
    ess.MoE.nu[i, ] <- coda::effectiveSize(as.mcmc(fit$MoEfitBayes$nu[-c(1:burnin),BayesOrder]))
    
    ## extract the first Sigma's entry from each of the 3 components
    Sigmas <- cbind(unlist(lapply(fit$MoEfitBayes$Sigma[-c(1:burnin)], function(xx) xx[1, 1, BayesOrder[1]])), 
                    unlist(lapply(fit$MoEfitBayes$Sigma[-c(1:burnin)], function(xx) xx[1, 1, BayesOrder[2]])), 
                    unlist(lapply(fit$MoEfitBayes$Sigma[-c(1:burnin)], function(xx) xx[1, 1, BayesOrder[3]])) )
    ess.MoE.Sigma[i, ] <- coda::effectiveSize(as.mcmc(Sigmas))
    
    if( dim(fit$MoEfitBayes$Beta_samples)[2] > 1 ) {
      # extract the first coefficients from two components; always the last dimension's value 0
      Betas <- fit$MoEfitBayes$Beta_samples[-c(1:burnin), 1, 1:2] 
      ess.MoE.beta[i, ] <- coda::effectiveSize(as.mcmc(Betas))
    }
    
  }
  
  ### calculate Mean and Standard Deviations
  value <- paste(c(
    paste(paste0(round(colMeans(ess.MM.nu), 0), " (", round(apply(ess.MM.nu, 2, sd), 0), ")")[1:2], collapse = " & "),
    paste(paste0(round(colMeans(ess.MM.Sigma), 0), " (", round(apply(ess.MM.Sigma, 2, sd), 0), ")")[1:2], collapse = " & "),
    paste(paste0(round(colMeans(ess.MM.pi), 0), " (", round(apply(ess.MM.pi, 2, sd), 0), ")")[1:2], collapse = " & ")
  ), collapse = " && ")
  ESS.MM.sd <- rbind(ESS.MM.sd, data.frame(p = p, n = paste0("n=",n), ESS = value))
  
  value <- paste(c(
    paste(paste0(round(colMeans(ess.MoE.nu), 0), " (", round(apply(ess.MoE.nu, 2, sd), 0), ")")[1:2], collapse = " & "),
    paste(paste0(round(colMeans(ess.MoE.Sigma), 0), " (", round(apply(ess.MoE.Sigma, 2, sd), 0), ")")[1:2], collapse = " & ")
  ), collapse = " && ")
  ESS.MoE.sd <- rbind(ESS.MoE.sd, data.frame(p = p, n = paste0("n=",n), ESS = value))

}


ESS_trueMM <- list(ESS.MM.sd = ESS.MM.sd, ESS.MoE.sd = ESS.MoE.sd)
ESS_trueMM
