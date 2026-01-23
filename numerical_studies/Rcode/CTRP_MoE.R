rm(list = ls())

# user-defined directory to save simulation results as data files
file_directory <- "../moewishart/"
filename_save <- file_directory

# load preprocessed CTRP drug sensitivity's covariance observations
load(paste0(filename_save, "CTRP_viability_S_list.RData"))


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
## extract 1st component score as covariate of MoE model
x <- cbind(as.numeric(status)-1, pca_res$x[, 1])


drug_class <- NULL
library(moewishart)

# run Bayesian mixture model, K = 2
my_K_intial <- 2
set.seed(123)
fitBayes <- moewishartX(
  S_list, X = cbind(1, x),
  K = my_K_intial, 
  nu_prior_a = 2, nu_prior_b = 0.2,
  mh_sigma = c(0.05, 0.1), mh_beta = 0.1,
  niter = 50000, burnin = 10000, thin = 1, verbose = TRUE
)
# check elpd-loo for tuning hyperparameter
nIter <- 10000; burnin <- 1000; thin <- 1
(a <- loo::loo(fitBayes$loglik_individual[((1 + burnin):nIter) / thin, ]))
fitBayes$elpd <- a$estimates["elpd_loo",]
save(fitBayes, file = paste0(filename_save, "MoEfitBayesK", my_K_intial, ".RData"))

z_hat <- round(colMeans(fitBayes$z[-c(1:10000),]), 0)
drug_class <- data.frame(drug_class, Bayes_MoE_K2 = z_hat)


# run Bayesian mixture model, K = 3
my_K_intial <- 3
set.seed(123)
fitBayes <- moewishartX(
  S_list, X = cbind(1, x),
  K = my_K_intial, 
  nu_prior_a = 2, nu_prior_b = 0.2,
  mh_sigma = c(0.07, 0.07, 0.1), mh_beta = 0.08,
  niter = 50000, burnin = 10000, thin = 1, verbose = TRUE
)

# check elpd-loo for tuning hyperparameter
nIter <- 10000; burnin <- 1000; thin <- 1
(a <- loo::loo(fitBayes$loglik_individual[((1 + burnin):nIter) / thin, ]))
fitBayes$elpd <- a$estimates["elpd_loo",]
save(fitBayes, file = paste0(filename_save, "MoEfitBayesK", my_K_intial, ".RData"))

z_hat <- round(colMeans(fitBayes$z[-c(1:10000),]), 0)
drug_class <- cbind(drug_class, Bayes_MoE_K3 = z_hat)


# run Bayesian mixture model, K = 4
my_K_intial <- 4
set.seed(123)
fitBayes <- moewishartX(
  S_list, X = cbind(1, x),
  K = my_K_intial, 
  nu_prior_a = 2, nu_prior_b = 0.2,
  mh_sigma = c(0.07, 0.12, 0.2, 0.07), 
  mh_beta = c(0.08, 0.1, 0.17),
  niter = 50000, burnin = 10000, thin = 1, verbose = TRUE
)
# check elpd-loo for tuning hyperparameter
nIter <- 10000; burnin <- 1000; thin <- 1
(a <- loo::loo(fitBayes$loglik_individual[((1 + burnin):nIter) / thin, ]))
fitBayes$elpd <- a$estimates["elpd_loo",]
save(fitBayes, file = paste0(filename_save, "MoEfitBayesK", my_K_intial, ".RData"))

z_hat <- round(colMeans(fitBayes$z[-c(1:10000),]), 0)
drug_class <- cbind(drug_class, Bayes_MoE_K4 = z_hat)

# run Bayesian mixture model, K = 5
my_K_intial <- 5
set.seed(123)
fitBayes <- moewishartX(
  S_list, X = cbind(1, x),
  K = my_K_intial, 
  nu_prior_a = 2, nu_prior_b = 0.2,
  mh_sigma = c(0.08, 0.15, 0.2, 0.08, 0.13), 
  mh_beta = c(0.08, 0.15, 0.15, 0.08), 
  niter = 50000, burnin = 10000, thin = 1, verbose = TRUE
)
# check elpd-loo for tuning hyperparameter
nIter <- 10000; burnin <- 1000; thin <- 1
(a <- loo::loo(fitBayes$loglik_individual[((1 + burnin):nIter) / thin, ]))
fitBayes$elpd <- a$estimates["elpd_loo",]
save(fitBayes, file = paste0(filename_save, "MoEfitBayesK", my_K_intial, ".RData"))

z_hat <- round(colMeans(fitBayes$z[-c(1:10000),]), 0)
drug_class <- cbind(drug_class, Bayes_MoE_K5 = z_hat)


# run Bayesian mixture model, K = 6
my_K_intial <- 6
set.seed(123)
fitBayes <- moewishartX(
  S_list, X = cbind(1, x),
  K = my_K_intial, 
  nu_prior_a = 2, nu_prior_b = 0.2,
  mh_sigma = c(0.08, 2, 0.18, 0.3, 0.07, 0.1),  
  mh_beta = c(0.08, 2, 0.15, 0.18, 0.1), 
  niter = 50000, burnin = 10000, thin = 1, verbose = TRUE
)
# check elpd-loo for tuning hyperparameter
nIter <- 10000; burnin <- 1000; thin <- 1
(a <- loo::loo(fitBayes$loglik_individual[((1 + burnin):nIter) / thin, ]))
fitBayes$elpd <- a$estimates["elpd_loo",]
save(fitBayes, file = paste0(filename_save, "MoEfitBayesK", my_K_intial, ".RData"))

z_hat <- round(colMeans(fitBayes$z[-c(1:10000),]), 0)
drug_class <- cbind(drug_class, Bayes_MoE_K6 = z_hat)


# run Bayesian mixture model, K = 7
my_K_intial <- 7
set.seed(123)
fitBayes <- moewishartX(
  S_list, X = cbind(1, x),
  K = my_K_intial, 
  nu_prior_a = 2, nu_prior_b = 0.2,
  mh_sigma = c(0.08, 2.2, 0.18, 0.2, 0.11, 0.09, 3),  
  mh_beta = c(0.08, 2, 0.18, 0.15, 0.12, 0.1), 
  niter = 50000, burnin = 10000, thin = 1, verbose = TRUE
)
# check elpd-loo for tuning hyperparameter
nIter <- 10000; burnin <- 1000; thin <- 1
(a <- loo::loo(fitBayes$loglik_individual[((1 + burnin):nIter) / thin, ]))
fitBayes$elpd <- a$estimates["elpd_loo",]
save(fitBayes, file = paste0(filename_save, "MoEfitBayesK", my_K_intial, ".RData"))

z_hat <- round(colMeans(fitBayes$z[-c(1:10000),]), 0)
drug_class <- cbind(drug_class, Bayes_MoE_K7 = z_hat)


# run Bayesian mixture model, K = 8
my_K_intial <- 8
set.seed(123)
fitBayes <- moewishartX(
  S_list, X = cbind(1, x),
  K = my_K_intial, #init_pi = init_pi, 
  nu_prior_a = 2, nu_prior_b = 0.2,
  mh_sigma = c(2.5, 3, 0.08, 0.2, 0.07, 0.15, 0.15, 3), 
  mh_beta = c(2, 1.8, 0.08, 0.2, 0.1, 0.18, 0.1), 
  niter = 50000, burnin = 10000, thin = 1, verbose = TRUE
)
# check elpd-loo for tuning hyperparameter
nIter <- 10000; burnin <- 1000; thin <- 1
(a <- loo::loo(fitBayes$loglik_individual[((1 + burnin):nIter) / thin, ]))
fitBayes$elpd <- a$estimates["elpd_loo",]
save(fitBayes, file = paste0(filename_save, "MoEfitBayesK", my_K_intial, ".RData"))

z_hat <- round(colMeans(fitBayes$z[-c(1:10000),]), 0)
drug_class <- cbind(drug_class, Bayes_MoE_K8 = z_hat)
metaData_selected_classes <- cbind(metaData_selected, status=status, drug_class)
write.csv(metaData_selected_classes, file = paste0(filename_save, "drug_class_MoE.csv"))

nIter <- 50000; burnin <- 0 # trace plot also shows burn-in period
datSim <- data.frame(K = NULL, iter = NULL, loglik = NULL, component = NULL, value = NULL)

## MCMC diagnosis:
for (K in 2:8) {

  load(paste0(filename_save, "MoEfitBayesK", K, ".RData"))
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

  betas <- fit$Beta_samples[, 2, -K] # remove intercepts and (last) reference group
  datSim <- rbind(datSim, data.frame(K = paste0("K=",K), iter = 1:burnin, component = NA, estimator = "beta", value = NA))
  datSim <- rbind(datSim, data.frame(K = paste0("K=",K), iter = rep(burnin+1:(nIter-burnin), K-1), component = rep(1:(K-1), each=nIter-burnin), estimator = "beta", value = as.vector(betas)))

}
save(datSim, file = paste0(filename_save, "CTRP_datSimMoE.RData"))

# summarize posterior mean and 95% credible interval
dat_Estimates <- data.frame(K = NULL, elpd = NULL, nu1 = NULL, nu2 = NULL, sigma1 = NULL, sigma1 = NULL, beta1 = NULL, beta2 = NULL)
nIter <- 50000; burnin <- 10000
digits0 <- 2
for (K in 2:8) {

  load(paste0(filename_save, "MoEfitBayesK", K, ".RData"))
  fit <- fitBayes

  nu.mcmc <- fit$nu[-c(1:burnin), ]
  nu12 <- sapply(1:2, function(xx) {
      tmp <- nu.mcmc[,xx]
      m <- round(mean(tmp), digits = digits0)
      lb <- round(quantile(tmp, 0.025), digits = digits0)
      ub <- round(quantile(tmp, 0.975), digits = digits0)
      paste0(m, " (", lb, ",", ub, ")")
    })

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

  beta.mcmc <- fit$Beta_samples[-c(1:burnin), 2, 1:2]
  beta12 <- sapply(1:2, function(xx) {
      tmp <- beta.mcmc[,xx]
      m <- round(mean(tmp), digits = digits0)
      lb <- round(quantile(tmp, 0.025), digits = digits0)
      ub <- round(quantile(tmp, 0.975), digits = digits0)
      paste0(m, " (", lb, ",", ub, ")")
  })

  elpd <- round(fit$elpd, digits = digits0)
  elpd <- paste0(elpd[1], " (", elpd[2], ")")

  dat_Estimates <- rbind(dat_Estimates, 
    data.frame(K = K, elpd = elpd, nu1 = nu12[1], nu2 = nu12[2], sigma1 = sigma12[1], sigma1 = sigma12[2], beta1 = beta12[1], beta2 = beta12[2])
  )

}
write.csv(dat_Estimates, file = paste0(filename_save, "drug_MoE_estimates.csv"))
paste0(dat_Estimates$elpd, collapse = " & ")


load(paste0(filename_save, "CTRP_datSimMoE.RData"))

datSim$K <- factor(datSim$K, levels = paste0("K=", 2:8))
datSim$estimator <- factor(datSim$estimator, 
  levels = c("loglik", "beta", "nu", "log|Sigma|"),
  labels = c("Log-likelihood", "beta", "nu", "log*'|'*Sigma*'|'"))
datSim$component <- factor(datSim$component)
library(ggplot2)
library(ggh4x)

theme_set(theme_bw())
gg <- ggplot(datSim, aes(y=value, x=iter, group = component, color = component)) + 
  geom_line(alpha=.8) + 
  ggtitle("MoE model") + 
  facet_grid(estimator ~ K, scales = "free", labeller = labeller(estimator = label_parsed, K = label_value)) +
  scale_color_discrete(breaks = as.character(1:K)) + 
  guides(color = guide_legend(title = "Component")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) #+
  # theme_bw()

gg <- gg + xlab("Iteration") + ylab("Value")

pdf("moewishart_CTRP_MoE.pdf", height = 6.5, width = 9)
gg
dev.off()



# load CTRP drug sensitivity's AUC data
# NOTE: mannually copy the cluster indices from four methods into the CSV file `metaData_selected.csv` and save it as `metaData_selected2.csv`
aucData <- read.csv(paste0(filename_save, "v20.data.curves_post_qc.txt"), header=T, sep="\t", fill=T)

load(paste0(filename_save, "CTRP_viability.RData"))
length(unique(viability_selected$master_cpd_id))
auc <- aucData[aucData$master_cpd_id %in% unique(viability_selected$master_cpd_id), ]
auc <- auc[auc$experiment_id %in% unique(viability_selected$experiment_id), ]

auc2 <- auc[, c(1,17,16)]
auc3 <- reshape(auc2, v.names="area_under_curve", timevar="master_cpd_id", idvar="experiment_id", direction="wide")
auc3[1:5,1:4]
length(unique(auc3$experiment_id))
sum(is.na(auc3))
colSums(is.na(auc3))


#===============
# select nonmissing AUC data
#===============
y0 <- auc3[, -1]
m0 <- dim(y0)[2]-1
eps <- 0.05
# r1.na is better to be not smaller than r2.na
r1.na <- 0.3
r2.na <- 0.2
k <- 1
while(sum(is.na(y0[,2:(1+m0)]))>0){
  r1.na <- r1.na - eps/k
  r2.na <- r1.na - eps/k
  k <- k + 1
  ## select drugs with <30% (decreasing with k) missing data overall cell lines
  na.y <- apply(y0[,2:(1+m0)], 2, function(xx) sum(is.na(xx))/length(xx))
  while(sum(na.y<r1.na)<m0){
    y0 <- y0[,-c(1+which(na.y>=r1.na))]
    m0 <- sum(na.y<r1.na)
    na.y <- apply(y0[,2:(1+m0)], 2, function(xx) sum(is.na(xx))/length(xx))
  }
  
  ## select cell lines with treatment of at least 80% (increasing with k) drugs
  na.y0 <- apply(y0[,2:(1+m0)], 1, function(xx) sum(is.na(xx))/length(xx))
  while(sum(na.y0<r2.na)<(dim(y0)[1])){
    y0 <- y0[na.y0<r2.na,]
    na.y0 <- apply(y0[,2:(1+m0)], 1, function(xx) sum(is.na(xx))/length(xx))
  }
  num.na <- sum(is.na(y0[,2:(1+m0)]))
  message("#{NA}=", num.na, "\n", "r1.na =", r1.na, ", r2.na =", r2.na, "\n")
}

dim(y0)
master_cpd_id172 <- substr(colnames(y0), 18, nchar(colnames(y0)))
y1 <- y0
colnames(y1) <- master_cpd_id172

# load meta data together with four methods' grouping results
metaData_selected_grouped <- read.csv(paste0(filename_save, "metaData_selected2.csv"), header=T, sep=";", fill=T)
metaData_selected_grouped172 <- metaData_selected_grouped[metaData_selected_grouped$master_cpd_id %in% colnames(y1), ]

y2 <- y1[, match(metaData_selected_grouped172$master_cpd_id, colnames(y1))]
sum(colnames(y2) == metaData_selected_grouped172$master_cpd_id)


metaData_grouped172 <- metaData_selected_grouped172[, c(2, 3, 13:18)]
length(unique(metaData_grouped172$EM_K8))

names(metaData_grouped172)[names(metaData_grouped172) == "combined_group"] <- "MoA"
metaData_grouped172$MoA[metaData_grouped172$MoA == "Miscellaneous or screening hit/natural product"] <- "Miscellaneous/screening hit"
metaData_grouped172$MoA <- factor(metaData_grouped172$MoA)

## plot heatmap with annotations
library(ComplexHeatmap)
anno <- metaData_grouped172
anno$status <- factor(anno$status, levels = c("approved_or_clinical", "experimental"),
  labels = c("FDA approved or clinical", "experimental"))
names(anno)[3:7] <- c("Status", "EM", "EM-MoE", "Bayes", "Bayes-MoE")
rownames(anno) <- anno$master_cpd_id
anno <- anno[, -c(1:2)]

# Match the labels between different methods
EM_k1 <- which(anno$`EM` == 1)
EM_k2 <- which(anno$`EM` == 2)
anno[EM_k1,"EM"] <- 2
anno[EM_k2,"EM"] <- 1
EM_k2 <- which(anno$`EM` == 2)
EM_k3 <- which(anno$`EM` == 3)
anno[EM_k2,"EM"] <- 3
anno[EM_k3,"EM"] <- 2
EM_k5 <- which(anno$`EM` == 5)
EM_k8 <- which(anno$`EM` == 8)
anno[EM_k5,"EM"] <- 8
anno[EM_k8,"EM"] <- 5
EM_k6 <- which(anno$`EM` == 6)
EM_k8 <- which(anno$`EM` == 8)
anno[EM_k6,"EM"] <- 8
anno[EM_k8,"EM"] <- 6

EM_MoE_k2 <- which(anno$`EM-MoE` == 2)
EM_MoE_k5 <- which(anno$`EM-MoE` == 5)
anno[EM_MoE_k2,"EM-MoE"] <- 5
anno[EM_MoE_k5,"EM-MoE"] <- 2
EM_MoE_k4 <- which(anno$`EM-MoE` == 4)
EM_MoE_k6 <- which(anno$`EM-MoE` == 6)
anno[EM_MoE_k4,"EM-MoE"] <- 6
anno[EM_MoE_k6,"EM-MoE"] <- 4
EM_MoE_k5 <- which(anno$`EM-MoE` == 5)
EM_MoE_k8 <- which(anno$`EM-MoE` == 8)
anno[EM_MoE_k5,"EM-MoE"] <- 8
anno[EM_MoE_k8,"EM-MoE"] <- 5

y3 <- scale(y2)
rownames(y3) <- paste0("cell_line", 1:nrow(y3))

# Define distinct color sets
status_colors = c("FDA approved or clinical" = "#E41A1C", "experimental" = "#377EB8")

combined_levels <- levels(as.factor(anno$MoA)) 
combined_colors <- setNames(ggsci::pal_igv()(11), combined_levels)
names(combined_colors) <- combined_levels

# Define Shared Colors 
shared_levels <- 1:10 
shared_colors <- structure(RColorBrewer::brewer.pal(10, "Paired"), names = shared_levels)

ordered_cols <- c("Status", "MoA", "Bayes", "Bayes-MoE", "EM", "EM-MoE")
anno <- anno[, ordered_cols]
ann_cols <- list(Status = status_colors, MoA = combined_colors)
for(name in setdiff(colnames(anno), c("Status", "MoA"))) {
    ann_cols[[name]] = shared_colors
}

# Define heatmap
col_ann <- HeatmapAnnotation(df = anno, col = ann_cols, show_legend = FALSE)
col_fun = colorRamp2(c(min(y3), 0, max(y3)), c("blue", "white", "red"))

ht <- Heatmap(y3, 
              name = "Z-Score", 
              col = col_fun,
              top_annotation = col_ann,
              show_row_names = FALSE, 
              show_column_names = FALSE,
              cluster_rows = FALSE,      
              show_heatmap_legend = FALSE, 
              column_title = "Compounds", 
              column_title_side = "bottom",
              row_title = "Cell Lines")

# Legends
lgd_zscore = Legend(title = "Standardized\nAUC", col_fun = col_fun, 
                    title_gp = gpar(fontsize = 8, fontface = "bold"),
                    legend_height = unit(2.5, "cm"))
lgd_status = Legend(title = "Status", at = names(status_colors), 
                    legend_gp = gpar(fill = status_colors),
                    grid_height = unit(3, "mm"), grid_width = unit(3, "mm"),
                    title_gp = gpar(fontsize = 11, fontface = "bold"),
                    labels_gp = gpar(fontsize = 10))
lgd_moa = Legend(title = "MoA", at = names(combined_colors), 
                      legend_gp = gpar(fill = combined_colors),
                      ncol = 1, # Split into two columns of 7
                      grid_height = unit(3, "mm"), grid_width = unit(3, "mm"),
                      title_gp = gpar(fontsize = 9, fontface = "bold"),
                      labels_gp = gpar(fontsize = 8))
lgd_shared = Legend(title = "Component index", at = names(shared_colors), 
                    legend_gp = gpar(fill = shared_colors),
                    grid_height = unit(3, "mm"), grid_width = unit(3, "mm"),
                    title_gp = gpar(fontsize = 11, fontface = "bold"),
                    labels_gp = gpar(fontsize = 10),
                    ncol = 1) # Set to 2 if too tall


pdf("moewishart_CTRP_heatmap.pdf",width=6, height=6)

# Draw heatmap with wide right padding for 3 columns of legends
draw(ht, padding = unit(c(2, 2, 2, 53), "mm")) 
draw(lgd_zscore, x = unit(0.52, "npc"), y = unit(0.55, "npc"), just = "left")
draw(lgd_status, x = unit(0.7, "npc"), y = unit(0.83, "npc"), just = "left")
draw(lgd_moa, x = unit(0.52, "npc"), y = unit(0.23, "npc"), just = "left")
draw(lgd_shared, x = unit(0.7, "npc"), y = unit(0.59, "npc"), just = "left")

dev.off()