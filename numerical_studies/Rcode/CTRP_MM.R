
rm(list = ls())

# user-defined directory to save simulation results as data files
file_directory <- "../moewishart/"
filename_save <- file_directory

# load preprocessed CTRP drug sensitivity's covariance observations
load(paste0(filename_save, "CTRP_viability_S_list.RData"))

library(moewishart)

# run Bayesian mixture model, K = 2
my_K_intial <- 2
set.seed(123)
fitBayes <- mixturewishart(S_list,
  K = my_K_intial, 
  nu_prior_a = 6, nu_prior_b = 0.5,
  mh_sigma = 0.06,
  niter = 50000, burnin = 10000, thin = 1, verbose = TRUE
)
# check elpd-loo for tuning hyperparameter
nIter <- 50000; burnin <- 10000; thin <- 1
(a <- loo::loo(fitBayes$loglik_individual[((1 + burnin):nIter) / thin, ]))
fitBayes$elpd <- a$estimates["elpd_loo",]
save(fitBayes, file = paste0(filename_save, "fitBayesK", my_K_intial, ".RData"))

z_hat <- round(colMeans(fitBayes$z[-c(1:10000),]), 0)
drug_class <- data.frame(Drugs = names(S_list), BayesMM_K2 = z_hat)


# run Bayesian mixture model, K = 3
my_K_intial <- 3
set.seed(123)
fitBayes <- mixturewishart(S_list,
  K = my_K_intial, 
  nu_prior_a = 2, nu_prior_b = 0.2,
  mh_sigma = 0.08,
  niter = 50000, burnin = 10000, thin = 1, verbose = TRUE
)
# check elpd-loo for tuning hyperparameter
nIter <- 50000; burnin <- 10000; thin <- 1
(a <- loo::loo(fitBayes$loglik_individual[((1 + burnin):nIter) / thin, ]))
fitBayes$elpd <- a$estimates["elpd_loo",]
save(fitBayes, file = paste0(filename_save, "fitBayesK", my_K_intial, ".RData"))

z_hat <- round(colMeans(fitBayes$z[-c(1:10000),]), 0)
drug_class <- cbind(drug_class, BayesMM_K3 = z_hat)


# run Bayesian mixture model, K = 4
my_K_intial <- 4
set.seed(123)
fitBayes <- mixturewishart(S_list,
  K = my_K_intial, 
  nu_prior_a = 2, nu_prior_b = 0.3,
  mh_sigma = c(0.08, 0.12, 0.2, 0.08),
  niter = 50000, burnin = 10000, thin = 1, verbose = TRUE
)
# check elpd-loo for tuning hyperparameter
nIter <- 50000; burnin <- 10000; thin <- 1
(a <- loo::loo(fitBayes$loglik_individual[((1 + burnin):nIter) / thin, ]))
fitBayes$elpd <- a$estimates["elpd_loo",]
save(fitBayes, file = paste0(filename_save, "fitBayesK", my_K_intial, ".RData"))

z_hat <- round(colMeans(fitBayes$z[-c(1:10000),]), 0)
drug_class <- cbind(drug_class, BayesMM_K4 = z_hat)


# run Bayesian mixture model, K = 5
my_K_intial <- 5
set.seed(123)
fitBayes <- mixturewishart(S_list,
  K = my_K_intial, 
  nu_prior_a = 2, nu_prior_b = 0.2,
  mh_sigma = c(0.08, 0.18, 0.2, 0.08, 0.1),
  niter = 50000, burnin = 10000, thin = 1, verbose = TRUE
)
# check elpd-loo for tuning hyperparameter
nIter <- 50000; burnin <- 10000; thin <- 1
(a <- loo::loo(fitBayes$loglik_individual[((1 + burnin):nIter) / thin, ]))
fitBayes$elpd <- a$estimates["elpd_loo",]
save(fitBayes, file = paste0(filename_save, "fitBayesK", my_K_intial, ".RData"))

z_hat <- round(colMeans(fitBayes$z[-c(1:10000),]), 0)
drug_class <- cbind(drug_class, BayesMM_K5 = z_hat)


# run Bayesian mixture model, K = 6
my_K_intial <- 6
set.seed(123)
fitBayes <- mixturewishart(S_list,
  K = my_K_intial, 
  nu_prior_a = 4, nu_prior_b = 0.5,
  mh_sigma = c(0.07, 2.2, 0.2, 0.15, 0.1, 0.2),
  niter = 50000, burnin = 10000, thin = 1, verbose = TRUE
)
# check elpd-loo for tuning hyperparameter
nIter <- 50000; burnin <- 10000; thin <- 1
(a <- loo::loo(fitBayes$loglik_individual[((1 + burnin):nIter) / thin, ]))
fitBayes$elpd <- a$estimates["elpd_loo",]
save(fitBayes, file = paste0(filename_save, "fitBayesK", my_K_intial, ".RData"))

z_hat <- round(colMeans(fitBayes$z[-c(1:10000),]), 0)
drug_class <- cbind(drug_class, BayesMM_K6 = z_hat)


# run Bayesian mixture model, K = 7
my_K_intial <- 7
set.seed(123)
fitBayes <- mixturewishart(S_list,
  K = my_K_intial, 
  nu_prior_a = 6, nu_prior_b = 0.5,
  mh_sigma = c(0.2, 2, 1.6, 0.08, 0.1, 0.08, 0.2),
  niter = 50000, burnin = 10000, thin = 1, verbose = TRUE
)
# check elpd-loo for tuning hyperparameter
nIter <- 50000; burnin <- 10000; thin <- 1
(a <- loo::loo(fitBayes$loglik_individual[((1 + burnin):nIter) / thin, ]))
fitBayes$elpd <- a$estimates["elpd_loo",]
save(fitBayes, file = paste0(filename_save, "fitBayesK", my_K_intial, ".RData"))

z_hat <- round(colMeans(fitBayes$z[-c(1:10000),]), 0)
drug_class <- cbind(drug_class, BayesMM_K7 = z_hat)


# run Bayesian mixture model, K = 8
my_K_intial <- 8
set.seed(123)
fitBayes <- mixturewishart(S_list,
  K = my_K_intial, 
  nu_prior_a = 2, nu_prior_b = 0.2,
  mh_sigma = c(2, 2.5, .05, 0.2, .08, 0.15, 0.15, 2.5),
  niter = 50000, burnin = 10000, thin = 1, verbose = TRUE
)
# check elpd-loo for tuning hyperparameter
nIter <- 50000; burnin <- 10000; thin <- 1
(a <- loo::loo(fitBayes$loglik_individual[((1 + burnin):nIter) / thin, ]))
fitBayes$elpd <- a$estimates["elpd_loo",]
save(fitBayes, file = paste0(filename_save, "fitBayesK", my_K_intial, ".RData"))

z_hat <- round(colMeans(fitBayes$z[-c(1:10000),]), 0)
drug_class <- cbind(drug_class, BayesMM_K8 = z_hat)
write.csv(drug_class, file = paste0(filename_save, "drug_class_MM.csv"))

nIter <- 50000; burnin <- 0 # trace plot also shows burn-in period
datSim <- data.frame(K = NULL, iter = NULL, loglik = NULL, component = NULL, value = NULL)

## MCMC diagnosis:
for (K in 2:8) {

  load(paste0(filename_save, "fitBayesK", K, ".RData"))
  fit <- fitBayes
  ## Mixture model's diagnosis
  datSim <- rbind(datSim, data.frame(K = paste0("K=",K), iter = 1:length(fit$loglik), component=NA, estimator = "loglik", value = fit$loglik))
  Sigma <- sapply(1:(nIter-burnin), function(xx){
    s0 <- fit$Sigma[[xx]]
    sapply(1:K, function(k){as.numeric(determinant(s0[,,k], logarithm = TRUE)$modulus)})
    # sapply(1:K, function(k){as.numeric(determinant(s0[[k]], logarithm = TRUE)$modulus)}) ## for Rcpp code
  })
  datSim <- rbind(datSim, data.frame(K = paste0("K=",K), iter = 1:burnin, component=NA, estimator = "log|Sigma|", value = NA))
  datSim <- rbind(datSim, data.frame(K = paste0("K=",K), iter = rep(burnin+1:(nIter-burnin), K), component = rep(1:K, each=nIter-burnin), estimator = "log|Sigma|", value = as.vector(t(Sigma))))

  nu <- as.vector(fit$nu)
  datSim <- rbind(datSim, data.frame(K = paste0("K=",K), iter = 1:burnin, component = NA, estimator = "nu", value = NA))
  datSim <- rbind(datSim, data.frame(K = paste0("K=",K), iter = rep(burnin+1:(nIter-burnin), K), component = rep(1:K, each=nIter-burnin), estimator = "nu", value = nu))

  pi <- as.vector(fit$pi)
  datSim <- rbind(datSim, data.frame(K = paste0("K=",K), iter = 1:burnin, component = NA, estimator = "pi", value = NA))
  datSim <- rbind(datSim, data.frame(K = paste0("K=",K), iter = rep(burnin+1:(nIter-burnin), K), component = rep(1:K, each=nIter-burnin), estimator = "pi", value = pi))

}

save(datSim, file = paste0(filename_save, "CTRP_datSimMM.RData"))


# summarize posterior mean and 95% credible interval
dat_Estimates <- data.frame(K = NULL, elpd = NULL, nu1 = NULL, nu2 = NULL, sigma1 = NULL, sigma2 = NULL, pi1 = NULL, pi2 = NULL)
nIter <- 50000; burnin <- 10000
digits0 <- 2
for (K in 2:8) {

  load(paste0(filename_save, "fitBayesK", K, ".RData"))
  fit <- fitBayes

  nu.mcmc <- fit$nu[-c(1:burnin), ]
  nu12 <- sapply(1:2, function(xx) {
      tmp <- nu.mcmc[,xx]
      m <- round(mean(tmp), digits = digits0)
      lb <- round(quantile(tmp, 0.025), digits = digits0)
      ub <- round(quantile(tmp, 0.975), digits = digits0)
      paste0(m, " (", lb, ",", ub, ")")
  })

  sigma.mcmc <- sapply(1:(nIter-burnin), function(xx){
    s0 <- fit$Sigma[[xx]]
    s0[1, 1, 1:2]
  })
  sigma.mcmc <- t(sigma.mcmc)
  sigma12 <- sapply(1:2, function(xx) {
      tmp <- sigma.mcmc[,xx]
      m <- round(mean(tmp), digits = digits0)
      lb <- round(quantile(tmp, 0.025), digits = digits0)
      ub <- round(quantile(tmp, 0.975), digits = digits0)
      paste0(m, " (", lb, ",", ub, ")")
  })

  pi.mcmc <- fit$pi[-c(1:burnin), ]
  pi12 <- sapply(1:2, function(xx) {
      tmp <- pi.mcmc[,xx]
      m <- round(mean(tmp), digits = digits0)
      lb <- round(quantile(tmp, 0.025), digits = digits0)
      ub <- round(quantile(tmp, 0.975), digits = digits0)
      paste0(m, " (", lb, ",", ub, ")")
  })

  elpd <- round(fit$elpd, digits = digits0)
  elpd <- paste0(elpd[1], " (", elpd[2], ")")

  dat_Estimates <- rbind(dat_Estimates, 
    data.frame(K = K, elpd = elpd, nu1 = nu12[1], nu2 = nu12[2], sigma1 = sigma12[1], sigma2 = sigma12[2], pi1 = pi12[1], pi2 = pi12[2])
  )

}
write.csv(dat_Estimates, file = paste0(filename_save, "drug_MM_estimates.csv"))
at_Estimates <- read.csv(paste0(filename_save, "drug_MM_estimates.csv"), header = T)
paste0(dat_Estimates$elpd, collapse = " & ")

load(paste0(filename_save, "CTRP_datSimMM.RData"))

datSim$K <- factor(datSim$K, levels = paste0("K=", 2:8))
datSim$estimator <- factor(datSim$estimator, 
  levels = c("loglik", "pi", "nu", "log|Sigma|"),
  labels = c("Log-likelihood", "pi", "nu", "log*'|'*Sigma*'|'"))
datSim$component <- factor(datSim$component)
library(ggplot2)
library(ggh4x)

theme_set(theme_bw())
gg <- ggplot(datSim, aes(y=value, x=iter, group = component, color = component)) + 
  geom_line(alpha=.8) + 
  ggtitle("Mixture model") + 
  facet_grid(estimator ~ K, scales = "free", labeller = labeller(estimator = label_parsed, K = label_value)) +
  scale_color_discrete(breaks = as.character(1:K)) + 
  guides(color = guide_legend(title = "Component")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) #+
  #theme_bw()

gg <- gg + xlab("Iteration") + ylab("Value")

pdf("moewishart_CTRP_MM.pdf", height = 6.5, width = 9)
gg
dev.off()

