<!-- badges: 
[![r-universe](https://ocbe-uio.r-universe.dev/badges/moewishart)](https://ocbe-uio.r-universe.dev/moewishart)
[![License](https://img.shields.io/badge/License-GPLv3-brightgreen.svg)](https://www.gnu.org/licenses/gpl-3.0)
-->
[![R-CMD-check](https://github.com/zhizuio/moewishart/workflows/R-CMD-check/badge.svg)](https://github.com/zhizuio/moewishart/actions)

## 1. Mixture-of-Experts Wishart models for covariance data

This R-package `moewishart` provides maximum likelihood estimation (MLE) and Bayesian estimation for the **Wishart mixture model** and the **Wishart mixture-of-experts** (**MoE-Wishart**) model. 
It implements four different inference algorithms for the two model:

- mixture model of Wishart distributions:
    - EM algorithm for finding the MLE;
    - Bayesian approach using Gibbs-within-MH sampling algorithm.
- Mixture-of-Expert model, in which the gating probabilities depend on covariates:
   - EM-MoE algorithm for finding the MLE;
   - Bayesian-MoE approach using a Gibbs-within-MH sampling algorithm.  


## 2. Installation

Install the latest development version from [GitHub](https://github.com/zhizuio/moewishart):


``` r
#library("devtools")
devtools::install_github("zhizuio/moewishart")
```

## 3. Example: data generating process


<br> 

Data simulation from a MoE-Wishart model:

- Sample size $n = 200$
- Dimension of the Wishart distribution $p = 2$
- Number of latent components $K = 3$
- $q=3$ gating covariates $\mathbf X = [x_{ij}] \in \mathbb R^{n\times q}$, $x_{ij}\sim\text{N}(0,1)$, $i=1,...,n$, $j=1,...,q$
- Fixed covariate effects $\boldsymbol\beta=[\boldsymbol\beta_1,...,\boldsymbol\beta_K] \in \mathbb R^{q\times K}$, with $\boldsymbol\beta_{K}=0$
- Probabilities of subpopulations $\boldsymbol\pi = [\pi_{ik}] \in \mathbb R^{n\times q}$, $\pi_{ik} = \exp(\mathbf X_i\boldsymbol\beta_k) / \sum_{l=1}^K\exp(\mathbf X_i\boldsymbol\beta_l)$, $k=1,...,K$
- Degrees of freedom $\boldsymbol\nu = (\nu_1,\nu_2,\nu_3) = (8, 12, 3)$
- Scale matrices of the Wishart distribution $\Sigma_1$, $\Sigma_2$, $\Sigma_3 \in \mathbb R^{p\times p}$
- Data $S_i \sim \pi_{i1}\text{Wishart}(\nu_1, \Sigma_1) + \pi_{i2}\text{Wishart}(\nu_2, \Sigma_2) + \pi_{i3}\text{Wishart}(\nu_3, \Sigma_3)$


``` r
library(moewishart)

n <- 200 # number of subjects
p <- 2 # dimension of covariance matrix
set.seed(123) # fix coefficients of underlying MoE model
Xq <- 3; K = 3
betas <- matrix(runif(Xq * K, -2, 2), nrow = Xq, ncol = K)
betas[, K] <- 0

# simulate data
dat <- simData(n, p,
  Xq = 3, K = 3, betas = betas,
  pis = c(0.35, 0.40, 0.25),
  nus = c(8, 12, 3)
)
```

## 4. Model fitting examples 

### 4.1 Bayesian MoE-Wishart model


``` r
# fit Bayesian MoE-Wishart model
set.seed(123)
fit <- moewishart(
  dat$S, X = cbind(1, dat$X), K = 3, 
  mh_sigma = c(0.2, 0.1, 0.2), # RW-MH variances (length K)
  mh_beta = c(0.3, 0.3), # RW-MH variances (length K-1)
  niter = 3000, burnin = 1000
)
```


### 4.2 Bayesian Wishart mixture model


``` r
# fit Bayesian Wishart mixture model
set.seed(123)
fit2 <- mixturewishart(
  dat$S, K = 3, 
  mh_sigma = c(0.2, 0.1, 0.2), # RW-MH variances
  niter = 3000, burnin = 1000
)
```


### 4.3 EM for MoE-Wishart


``` r
# fit MoE-Wishart model via EM alg.
set.seed(123)
fit3 <- moewishart(
  dat$S, X = cbind(1, dat$X), K = 3, 
  method = "em",
  niter = 3000
)
```


### 4.4 EM for Wishart mixture


``` r
# fit Wishart mixture model via EM alg.
set.seed(123)
fit4 <- mixturewishart(
  dat$S, K = 3, 
  method = "em",
  niter = 3000
)
```



