##===================
## Analyze data for paper by The Tien Mai & Zhi Zhao
##===================

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("PharmacoGx")
# 
# library(PharmacoGx) # PharmacoGx_3.12.2 
# availablePSets()
# # GDSC <- downloadPSet("GDSC_2020(v2-8.2)")
# CTRP <- downloadPSet("CTRPv2_2015")
# 
# a <- sensitivityRaw(CTRP)
# dim(a)


#===============
# ## load CTRP pharmacological data downloaded from https://ocg.cancer.gov/programs/ctd2/data-portal
#===============

# load drug viability
viability <- read.table("/Users/zhiz/Downloads/IPFStructPenalty/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.data.per_cpd_well.txt", header=T)
unique(viability$cpd_conc_umol)

# load compound names
compound <- read.table("/Users/zhiz/Downloads/IPFStructPenalty/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_compound.txt", fill=T)
compound <- compound[,1:2]
colnames(compound) <- compound[1,]

# merge drug viability and compound names
viability1 <- merge(viability, compound, by="master_cpd_id")


# remove viability from drug combinations; keep only monotherapy drug treatment
library(data.table)
viability1 <- viability1[!viability1$cpd_name %like% ":", ]

# filter drugs with common concentrations
concentration_common <- data.frame(drugID = viability1$master_cpd_id, 
                                   concentration = viability1$cpd_conc_umol)
concentration_common$drugID <- factor(concentration_common$drugID)
concentration_common$concentration <- factor(concentration_common$concentration)

# 496 drugs x 300 doses
tab <- table(concentration_common$drugID, concentration_common$concentration)
max(rowSums(tab > 0)) # one drug has at most 44 doses
min(rowSums(tab > 0)) # one drug has at least 16 doses

which.min(rowSums(tab > 0)) # drugID=1788 has the least 16 doses
sum(rowSums(tab > 0) == 16) # 475 drugs with each has 16 doses (may different over drugs)

tab16 <- tab[rowSums(tab > 0) == 16, ]
sum(tab16 == 0)

y <- data.matrix(tab16)
rownames(y) <- rownames(tab16)
y[y == 0] <- NA

#===============
# select nonmissing pharmacological data
#===============
y0 <- y
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

apply(y0, 2, function(xx) sum(is.na(xx)))

sum(is.na(y0[, -1]))
y <- y0[, -1]
# Complete 374 drugs x 15 doses
tab.full <- tab16[rownames(y), colnames(y)]
tab.full[1:4,1:3] 

# filter out small or large doses with little variation of viability, or remove every other doses
# keep ca. 4-8 doses

# keep every other doses
viability_selected <- viability1[viability1$cpd_conc_umol %in% colnames(tab.full)[seq(1,ncol(tab.full),by=3)] &
                                   viability1$master_cpd_id %in% rownames(tab.full), 
                                #c("assay_plate_barcode", "raw_value_log2", "cpd_conc_umol", "master_cpd_id")
                                ]
save(viability_selected, file = "CTRP_viability.RData")
#y <- reshape(viability_selected, v.names="raw_value_log2", timevar="assay_plate_barcode", idvar="master_cpd_id", direction="wide")

# compute covariance matrix of each drug based on replicates of each does
drugIDs <- unique(viability_selected$master_cpd_id)
drugNames <- unique(viability_selected$cpd_name)

n <- length(drugIDs)

# S <- list()
# for (i in 1:n) {
#   ## viability metric via dmso_zscore_log2: log2-transformed compound-well value z-scored using dmso_plate_avg_log2 and dmso_expt_std_log2
#   # dmso_plate_avg_log2:	per-plate mean of log2-transformed DMSO-well values
#   # dmso_expt_std_log2: per-experiment standard deviation of background-subtracted log2-transformed DMSO-well values
#   oneDrug <- viability_selected[viability_selected$master_cpd_id == drugIDs[i], 
#                                 c("assay_plate_barcode", "dmso_zscore_log2", "cpd_conc_umol")]
#   # reshape the data from long-format to wide-format
#   y <- reshape(oneDrug, v.names="dmso_zscore_log2", timevar="assay_plate_barcode", 
#                  idvar="cpd_conc_umol", direction="wide")
#   rownames(y) <- y[, 1]
#   y <- t(y[, -1])
#   # remove missing replicates
#   y.complete <- y[rowSums(is.na(y))==0, order(colnames(y))]
#   S[[i]] <- cov(y.complete)
#   names(S)[i] <- drugNames[i]
# }
# summary( unlist(lapply(S, det)) )

## If use raw viability (log2-scale directly)
S_list <- list()
for (i in 1:length(drugIDs)) {
  oneDrug <- viability_selected[viability_selected$master_cpd_id == drugIDs[i], c("assay_plate_barcode", "raw_value_log2", "cpd_conc_umol")]
  # reshape the data from long-format to wide-format
  y <- reshape(oneDrug, v.names="raw_value_log2", timevar="assay_plate_barcode", idvar="cpd_conc_umol", direction="wide")
  rownames(y) <- y[, 1]
  y <- t(y[, -1])
  # remove missing replicates
  y.complete <- y[rowSums(is.na(y))==0, order(colnames(y))]
  S_list[[i]] <- cov(y.complete)
  names(S_list)[i] <- drugNames[i]
}
summary( unlist(lapply(S_list, det)) )


##==========================
## Modeling covariances by mixture model and MoE model via Bayesian and EM algorithms
##==========================

my_K_intial <- 3
init_pi <- c(0.33, 0.33, 0.34)

# run Bayesian mixture model
set.seed(123)
fitBayes <- moewishart(S_list,
                       K = my_K_intial, #init_pi = init_pi, 
                       nu_prior_a = 4, nu_prior_b = 0.5,
                       mh_sigma = 0.07, #cpp = TRUE,
                       niter = 10000, burnin = 1000, thin = 1, verbose = TRUE
                       #niter = 20000, burnin = 5000, thin = 1, verbose = TRUE
)
# check elpd-loo for tuning hyperparameter
nIter <- 10000; burnin <- 1000; thin <- 1
loo::loo(fit$loglik_individual[((1 + burnin):nIter) / thin, ])

# run mixture model with EM algorithm
set.seed(123)
fitEM <- moewishart(S_list,
                    method = "em", cpp = TRUE,
                    K = my_K_intial, init_pi = init_pi, 
                    niter = 10000, verbose = TRUE
)

# run Bayesian mixture model
set.seed(123)
MoEfitBayes <- moewishartX(
  S_list, X = matrix(rep(1, n), ncol = 1),
  K = my_K_intial, #init_pi = init_pi, 
  nu_prior_a = 4, nu_prior_b = 0.5,
  mh_sigma = 0.07, mh_beta = 0.4,
  niter = 20000, burnin = 5000, thin = 1, verbose = TRUE
)

# run MoE model with EM algorithm
set.seed(123)
MoEfitEM <- moewishartX(
  S_list, X = matrix(rep(1, n), ncol = 1),
  method = "em",
  K = my_K_intial, #init_pi = init_pi, 
  niter = 10000, verbose = TRUE
)

# put together the results of all 4 methods
fit <- list(fitBayes = fitBayes, fitEM = fitEM, MoEfitBayes = MoEfitBayes, MoEfitEM = MoEfitEM)




