
rm(list = ls())

# user-defined directory to save simulation results as data files
file_directory <- "../moewishart/"
filename_save <- file_directory

# load preprocessed CTRP drug sensitivity's covariance observations
load(paste0(filename_save, "CTRP_viability_S_list.RData"))

library(moewishart)

# run EM mixture model
K <- 2:11 #c(2, 3, 5, 8, 10, 15, 20, 30)
fitEM <- list()
for(i in 1:length(K)) {
  set.seed(123)
  fitEM[[i]] <- mixturewishart(S_list,
    method = "em", #cpp = TRUE,
    K = K[i], 
    niter = 10000, verbose = TRUE,
    n_restarts = 20, restart_iters = 40, tol = 1e-8
  )
}

# load meta data of drugs
metaData <- read.csv(paste0(filename_save, "CTRP/v20.meta.per_compound.csv"), header=T, sep=";", fill=T)
metaData_selected <- metaData[metaData$master_cpd_id %in% names(S_list), ]
metaData_selected <- metaData_selected[match(names(S_list), metaData_selected$master_cpd_id), ]
status <- metaData_selected$cpd_status
status[status %in% c("FDA", "clinical")] <- "approved_or_clinical"
status[status %in% c("probe", "GE-active")] <- "experimental"
status <- factor(status)

# Fingerprint of SMILES
library(rcdk)
library(fingerprint)
mols <- parse.smiles(metaData_selected$cpd_smiles)
fps <- lapply(mols, get.fingerprint, type='circular')
fp.sim <- fingerprint::fp.sim.matrix(fps, method='tanimoto')
fp.dist <- 1 - fp.sim
# Perform PCA
pca_res <- prcomp(fp.dist, scale. = TRUE)
summ <- summary(pca_res)
head(summ$importance[2, ]) 
x <- cbind(as.numeric(status)-1, pca_res$x[, 1])


# run EM-MoE model
K <- 2:11 
fitEM_MoE <- list()
for(i in 1:length(K)) {
  set.seed(123)
  fitEM_MoE[[i]] <- moewishart(
    S_list, X = cbind(1, x),
    method = "em",
    K = K[i], 
    niter = 10000, verbose = TRUE
  )
}

fitEM_all <- list(fitEM = fitEM, fitEM_MoE = fitEM_MoE)
save(fitEM_all, file = paste0(filename_save, "CTRP_EM_all.RData"))

load(paste0(filename_save, "CTRP_EM_all.RData"))


# summarize EM estimates
p <- 5 # five drug doses
q <- 2 # one covariate per component
n <- length(S_list)
K <- 2:11
dat_Estimates <- data.frame(model = NULL, K = NULL, BIC = NULL, AIC = NULL, 
  nu1 = NULL, nu2 = NULL, sigma1 = NULL, sigma2 = NULL, pi_beta1 = NULL, pi_beta2 = NULL)
dat_IC <- data.frame(model = NULL, K = NULL, AIC = NULL, BIC = NULL, ICL = NULL)
drug_class <- NULL
drug_class <- data.frame(Drugs = names(S_list))

digits0 <- 2
for (i in 1:10) {

  # total number of estimated parameters
  kk <- K[i] * (p*(p+1)/2+1) + (K[i]-1) * (q+1)

  # results of mixture model via EM
  aic <- -2 * mean(fitEM[[i]]$loglik) + 2 * kk
  bic <- -2 * mean(fitEM[[i]]$loglik) + log(n) * kk
  nu12 <- round(fitEM[[i]]$nu, digits = digits0)
  sigma12 <- round(unlist(lapply(fitEM[[i]]$Sigma, function(xx) xx[1])), digits = digits0)
  pi12 <- round(fitEM[[i]]$pi, digits = digits0)

  dat_Estimates <- rbind(dat_Estimates, 
    data.frame(model = "EM", K = K[i], BIC = bic, AIC = aic, nu1 = nu12[1], nu2 = nu12[2], sigma1 = sigma12[1], sigma2 = sigma12[2], pi_beta1 = pi12[1], pi_beta2 = pi12[2])
  )

  r_ik <- fitEM[[i]]$tau
  icl <- bic - 2 * sum(r_ik * log(r_ik), na.rm = TRUE)
  dat_IC <- rbind(dat_IC, 
    data.frame(model = "EM", K = K[i], AIC = aic, BIC = bic, ICL = icl)
  )
  drug_class <- data.frame(drug_class, apply(r_ik, 1, which.max))
  names(drug_class)[ncol(drug_class)] <- paste0("EM_K", K[i])

  # results of MoE model via EM
  loglik <- fitEM_MoE[[i]]$loglik[length(fitEM_MoE[[i]]$loglik)]
  aic <- -2 * loglik + 2 * kk
  bic <- -2 * loglik + log(n) * kk
  nu12 <- round(fitEM_MoE[[i]]$nu, digits = digits0)
  sigma12 <- round(unlist(lapply(fitEM_MoE[[i]]$Sigma, function(xx) xx[1])), digits = digits0)
  beta12 <- round(fitEM_MoE[[i]]$Beta[1, ], digits = digits0)

  dat_Estimates <- rbind(dat_Estimates, 
    data.frame(model = "EM-MoE", K = K[i], BIC = bic, AIC = aic, nu1 = nu12[1], nu2 = nu12[2], sigma1 = sigma12[1], sigma2 = sigma12[2], pi_beta1 = beta12[1], pi_beta2 = beta12[2])
  )

  r_ik <- fitEM_MoE[[i]]$gamma
  icl <- bic - 2 * sum(r_ik * log(r_ik), na.rm = TRUE)
  dat_IC <- rbind(dat_IC, 
    data.frame(model = "EM-MoE", K = K[i], AIC = aic, BIC = bic, ICL = icl)
  )
  drug_class <- data.frame(drug_class, apply(r_ik, 1, which.max))
  names(drug_class)[ncol(drug_class)] <- paste0("EM_MoE_K", K[i])
}
dat_Estimates$model <- factor(dat_Estimates$model)
dat_Estimates[order(dat_Estimates$model), ]

dat_IC$model <- factor(dat_IC$model)
tmp <- sapply(which(dat_IC$model == "EM"), function(ii) t( dat_IC[ii, -c(1:2)] ) ); apply(tmp, 1, which.min) # EM
sapply(1:nrow(tmp), function(ii) paste0(round(tmp[ii,],0), collapse=" & "))
tmp <- sapply(which(dat_IC$model == "EM-MoE"), function(ii) t( dat_IC[ii, -c(1:2)] ) ); apply(tmp, 1, which.min) # EM
sapply(1:nrow(tmp), function(ii) paste0(round(tmp[ii,],0), collapse=" & "))
