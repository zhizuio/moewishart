rm(list = ls())

n <- 500 # subjects
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
  mh_sigma = 0.2, mh_beta = 0.4,
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

n <- 500
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
# n <- 500 # subjects
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

# source("moewishartX_errors.R")

# sapply(1:6, function(xx) as.vector(paste0(error.mean[xx, ], " (", error.sd[xx, ], ") ", collapse = "& ")))



# ########################
# ## summarize results of MoE models (model misspecification)
# ########################

# beta.true <- NULL
# source("moewishartX_errors.R")
# sapply(1:6, function(xx) as.vector(paste0(error.mean2[xx, ], " (", error.sd2[xx, ], ") ", collapse = "& ")))


