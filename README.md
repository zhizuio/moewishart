<!-- badges: 
[![r-universe](https://ocbe-uio.r-universe.dev/badges/moewishart)](https://ocbe-uio.r-universe.dev/moewishart)
[![License](https://img.shields.io/badge/License-GPLv3-brightgreen.svg)](https://www.gnu.org/licenses/gpl-3.0)
-->
[![R-CMD-check](https://github.com/zhizuio/moewishart/workflows/R-CMD-check/badge.svg)](https://github.com/zhizuio/moewishart/actions)

# moewishart: Mixture-of-Experts Wishart models for covariance data.

## The package contains 4 different algorithms for 2 model:
- mixture model of Wishart distributions:
    - EM algorithm for finding the MLE (maximum likelihood estimation);
    - Bayesian approach using Gibbs-within-MH sampling algorithm.
- Mixture-of-Expert model, in which the probability gating can depend on covariate(s):
   - EM-MoE algorithm for finding the MLE (maximum likelihood estimation);
   - Bayesian-MoE approach using Gibbs-within-MH sampling algorithm.  


## Installation

Install the latest development version from [GitHub](https://github.com/zhizuio/moewishart)

```r
#library("devtools")
devtools::install_github("zhizuio/moewishart", ref="main", auth_token = "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")
```

## Example
to be added, some usage example
