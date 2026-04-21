[![CRAN
status](https://www.r-pkg.org/badges/version/moewishart)](https://cran.r-project.org/package=moewishart)
[![r-universe](https://zhizuio.r-universe.dev/badges/moewishart)](https://zhizuio.r-universe.dev/moewishart)
[![R-CMD-check](https://github.com/zhizuio/moewishart/workflows/R-CMD-check/badge.svg)](https://github.com/zhizuio/moewishart/actions)
[![License](https://img.shields.io/badge/License-GPLv3-brightgreen.svg)](https://www.gnu.org/licenses/gpl-3.0)



This R-package `moewishart` provides maximum likelihood estimation (MLE) and Bayesian estimation for the **Wishart mixture model** and the **Wishart mixture-of-experts** (**MoE-Wishart**) model. 
It implements four different inference algorithms for the two model:

- mixture model of Wishart distributions:
    - EM algorithm for finding the MLE;
    - Bayesian approach using Gibbs-within-MH sampling algorithm.
- Mixture-of-Expert model, in which the gating probabilities depend on covariates:
   - EM-MoE algorithm for finding the MLE;
   - Bayesian-MoE approach using a Gibbs-within-MH sampling algorithm.  


## Installation

Install the latest released version from [CRAN](https://CRAN.R-project.org/package=moewishart):


``` r
install.packages("BayesSUR")
```

Install the latest development version from [GitHub](https://github.com/zhizuio/moewishart):


``` r
# library("devtools")
devtools::install_github("zhizuio/moewishart")
```

## Example


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


### 1. Working model: Bayesian MoE-Wishart model


``` r
library(moewishart)

n <- 200 # number of subjects
p <- 2 # dimension of covariance matrix
set.seed(123) # fix coefficients of underlying MoE model
Xq <- 3
K <- 3
betas <- matrix(runif(Xq * K, -2, 2), nrow = Xq, ncol = K)
betas[, K] <- 0

# simulate data
dat <- simData(n, p,
  Xq = 3, K = 3, betas = betas,
  pis = c(0.35, 0.40, 0.25),
  nus = c(8, 12, 3)
)

# fit Bayesian MoE-Wishart model
set.seed(123)
fit <- moewishart(
  dat$S,
  X = cbind(1, dat$X), K = 3,
  mh_sigma = c(0.2, 0.1, 0.2), # RW-MH variances (length K)
  mh_beta = c(0.3, 0.3), # RW-MH variances (length K-1)
  niter = 3000, burnin = 1000
)
```

<br> 

Posterior means for degrees of freedom (DoF) of Wishart distributions:


``` r
burnin <- 1000
nu_mcmc <- fit$nu[-c(1:burnin), ]
colMeans(nu_mcmc)
```

```
## [1]  8.574911 14.397351  3.310689
```

<br> 
True DoF:


``` r
dat$nu # true nu
```

```
## [1]  8 12  3
```

<br> 
Posterior means for scale matrices of Wishart distributions:


``` r
MoE_Sigma <- Reduce("+", fit$Sigma) / length(fit$Sigma)
MoE_Sigma
```

```
## , , 1
## 
##           [,1]      [,2]
## [1,] 0.5197160 0.2103881
## [2,] 0.2103881 0.7470847
## 
## , , 2
## 
##           [,1]      [,2]
## [1,] 1.7637949 0.5540576
## [2,] 0.5540576 1.3244947
## 
## , , 3
## 
##            [,1]       [,2]
## [1,]  4.1115070 -0.1267705
## [2,] -0.1267705  3.0385263
```

<br> 
Posterior means for gating coefficients:


``` r
beta_mcmc <- fit$Beta_samples[-c(1:burnin), , ]
apply(beta_mcmc, c(2, 3), mean)
```

```
##            [,1]        [,2] [,3]
## [1,] -0.3656861 -0.08024419    0
## [2,] -0.9526224  2.24956385    0
## [3,]  1.7609922  2.40287152    0
## [4,] -0.4953755 -2.56072719    0
```


### 2. Working model: Bayesian Wishart mixture model


``` r
# fit Bayesian Wishart mixture model
set.seed(123)
fit2 <- mixturewishart(
  dat$S,
  K = 3,
  mh_sigma = c(0.2, 0.1, 0.2), # RW-MH variances
  niter = 3000, burnin = 1000
)
```

<br> 
Posterior means for subpopulation probabilities:


``` r
colMeans(fit2$pi[-c(1:burnin), ])
```

```
## [1] 0.2690425 0.5088864 0.2220712
```

<br> 

Posterior means for DoF of Wishart distributions:


``` r
colMeans(fit2$nu[-c(1:burnin), ])
```

```
## [1]  7.986113 12.153338  3.284252
```


### 3. Working model: MoE-Wishart model via EM algorithm


``` r
# fit MoE-Wishart model via EM alg.
set.seed(123)
fit3 <- moewishart(
  dat$S,
  X = cbind(1, dat$X), K = 3,
  method = "em",
  niter = 3000
)
```

```
## Iter   1  loglik = -2079.322610
## Iter   2  loglik = -1998.694495
## Iter   3  loglik = -1985.659443
## Iter   4  loglik = -1947.223842
## Iter   5  loglik = -1899.666938
## Iter   6  loglik = -1878.233062
## Iter   7  loglik = -1861.702657
## Iter   8  loglik = -1851.548347
## Iter   9  loglik = -1847.342510
## Iter  10  loglik = -1845.497390
## Iter  11  loglik = -1844.693281
## Iter  12  loglik = -1844.360139
## Iter  13  loglik = -1844.220804
## Iter  14  loglik = -1844.160033
## Iter  15  loglik = -1844.132434
## Iter  16  loglik = -1844.119493
## Iter  17  loglik = -1844.113254
## Iter  18  loglik = -1844.110247
## Iter  19  loglik = -1844.108705
## Iter  20  loglik = -1844.107915
## Iter  21  loglik = -1844.107527
## Iter  22  loglik = -1844.107312
## Iter  23  loglik = -1844.107207
## Iter  24  loglik = -1844.107148
## Iter  25  loglik = -1844.107116
## Iter  26  loglik = -1844.107098
## Iter  27  loglik = -1844.107088
## Iter  28  loglik = -1844.107082
## Iter  29  loglik = -1844.107080
## Iter  30  loglik = -1844.107079
## Iter  31  loglik = -1844.107077
## Iter  32  loglik = -1844.107077
## Converged by loglik tolerance.
```

<br> 

EM estimates for DoF of Wishart distributions:


``` r
fit3$nu
```

```
## [1]  7.515417 13.987158  3.274665
```

<br> 
EM estimates for Wishart scale matrices:


``` r
fit3$Sigma
```

```
## [[1]]
##           [,1]      [,2]
## [1,] 0.5591113 0.2324429
## [2,] 0.2324429 0.8148737
## 
## [[2]]
##           [,1]      [,2]
## [1,] 1.7665723 0.5567668
## [2,] 0.5567668 1.3336367
## 
## [[3]]
##            [,1]       [,2]
## [1,]  4.3139885 -0.1886288
## [2,] -0.1886288  3.0983710
```

<br> 
EM estimates for gating coefficients:


``` r
fit3$Beta
```

```
##             comp1      comp2 comp3
## [1,] -0.006270492  0.1302039     0
## [2,] -0.798302303  2.0340525     0
## [3,]  1.598530103  2.2399293     0
## [4,] -0.510585695 -2.3465248     0
```


### 4. Working model: Wishart mixture model via EM algorithm


``` r
# fit Wishart mixture model via EM alg.
set.seed(123)
fit4 <- mixturewishart(
  dat$S,
  K = 3,
  method = "em",
  niter = 3000
)
```

```
## Running 3 initialization restarts...
##   -> Restart 1: Loglik = -2011.91
##   -> Restart 2: Loglik = -1989.16
##   -> Restart 3: Loglik = -1989.16
## Iter  10 | Loglik: -1953.0326 | Nu: 5, 7.26, 5
## Iter  20 | Loglik: -1952.7924 | Nu: 5.15, 7.05, 5.21
## Iter  30 | Loglik: -1930.6036 | Nu: 4.21, 7.13, 8.27
## Iter  40 | Loglik: -1902.9306 | Nu: 3.08, 7.78, 10.79
## Iter  50 | Loglik: -1902.6288 | Nu: 3, 7.76, 11.11
## Iter  60 | Loglik: -1902.6133 | Nu: 2.99, 7.71, 11.21
## Iter  70 | Loglik: -1902.6102 | Nu: 2.99, 7.69, 11.26
## Iter  80 | Loglik: -1902.6095 | Nu: 2.99, 7.68, 11.28
## Iter  90 | Loglik: -1902.6093 | Nu: 3, 7.68, 11.3
## Iter 100 | Loglik: -1902.6092 | Nu: 3, 7.68, 11.3
## Iter 110 | Loglik: -1902.6092 | Nu: 3, 7.68, 11.31
## Converged at iteration 116
```

<br> 
EM estimates for DoF of Wishart distributions:


``` r
fit4$nu
```

```
## [1]  2.995383  7.682819 11.309040
```

<br> 
EM estimate for Wishart scale matrices:


``` r
fit4$Sigma
```

```
## [[1]]
##           [,1]     [,2]
## [1,]  4.048859 -1.10103
## [2,] -1.101030  2.79641
## 
## [[2]]
##           [,1]      [,2]
## [1,] 0.5582012 0.2515063
## [2,] 0.2515063 0.8151752
## 
## [[3]]
##           [,1]      [,2]
## [1,] 2.0930529 0.6273737
## [2,] 0.6273737 1.5706730
```



