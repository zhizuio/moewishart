##========================================================================================
## This file contains code to produce Figure 2A. 
## 
## Input: Data files from R code 'simX_N200_p2.R', 'simX_N500_p2.R', 'simX_N1000_p2.R'
## Output: Boxplot
## 
## NOTE: To produce Figure 2B, 
##       just change the input data sets from "datSim_N200.RData" etc. to "datSim_X_N200.RData" etc.
##========================================================================================

rm(list = ls())

# specify folder directory including the results from R code 'simX_N200_p2.R', 'simX_N500_p2.R', 'simX_N1000_p2.R'
file_directory <- "../moewishart/"

# load 3 data sets from simulations with n=200, 500 and 1000
load(paste0(file_directory, "datSimX_N200.RData"))
load(paste0(file_directory, "datSimX_N500.RData"))
load(paste0(file_directory, "datSimX_N1000.RData"))

datSim <- rbind(datSimN200, datSimN500, datSimN1000)
datSim$n <- factor(datSim$n, levels = c("n=200", "n=500", "n=1000"))

datSim <- datSim[datSim$method %in% c("Bayes0", "EM", "Bayes0.X", "EM.X"), ]
datSim$method <- factor(datSim$method)
datSim$method <- factor(datSim$method, levels = c("Bayes0", "EM", "Bayes0.X", "EM.X"),
                        labels = c("Bayes", "EM", "Bayes-MoE", "EM-MoE"))

datSim <- datSim[datSim$estimator %in% c("nu-norm1", "Sigma-normF", "Beta-normF"), ]
datSim$estimator <- factor(datSim$estimator, levels = c("nu-norm1", "Sigma-normF", "Beta-normF"),
                           labels = c("frac(1,K)*'||'*hat(nu)-nu*'||'[1]", 
                                      "frac(1,K)*sum('||'*hat(Sigma)[k]-Sigma[k]*'||'[2]^2, k==1, K)", 
                                      "frac(1,K)*sum('||'*hat(beta)[k]-beta[k]*'||'[2]^2, k==1, K)"))

library(ggplot2)

theme_set(theme_bw())
my_colors <- c( "#d83034", "#ff9d3a", "#4ecb8d", "#008dff")
gg <- ggplot(datSim, aes(y=error, x=method)) + 
  geom_boxplot(aes(fill = method)) + 
  facet_grid(estimator ~ n, scales = "free", labeller = labeller(estimator = label_parsed, n = label_value)) +
  scale_fill_manual(values = my_colors) + 
  labs(x = "", y = "Error") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "none") 

gg
