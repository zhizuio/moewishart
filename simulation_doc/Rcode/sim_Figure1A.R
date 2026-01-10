##========================================================================================
## This file contains code to produce Figure 1A
## 
## Input: Data files from R code 'simN200_p2.R', 'simN500_p2.R', 'simN1000_p2.R'
## Output: Boxplot
##========================================================================================

rm(list = ls())

# specify folder directory including the results from R code 'simN200_p2.R', 'simN500_p2.R', 'simN1000_p2.R'
file_directory <- "../moewishart/"

# load 3 data sets from simulations with n=200, 500 and 1000
load(paste0(file_directory, "datSim_N200.RData"))
load(paste0(file_directory, "datSim_N500.RData"))
load(paste0(file_directory, "datSim_N1000.RData"))

datSim <- rbind(datSimN200, datSimN500, datSimN1000)
datSim$n <- factor(datSim$n, levels = c("n=200", "n=500", "n=1000"))

# exclude results from BMA and unnecessary metrics
datSim <- datSim[datSim$method %in% c("Bayes0", "EM", "Bayes0.X", "EM.X"), ]
datSim$method <- factor(datSim$method)
datSim$method <- factor(datSim$method, levels = c("Bayes0", "EM", "Bayes0.X", "EM.X"),
                        labels = c("Bayes", "EM", "Bayes-MoE", "EM-MoE"))

datSim <- datSim[datSim$estimator %in% c("nu-norm1", "pi-norm1", "Sigma-normF"), ]
# re-label estimators to show Greek letters
datSim$estimator <- factor(datSim$estimator, levels = c("pi-norm1", "nu-norm1", "Sigma-normF"),
                           labels = c("'||'*pi*'||'[1]", "'||'*nu*'||'[1]", "'||'*Sigma*'||'[2]"))

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
