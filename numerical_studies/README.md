# Numerical studies

This folder is for the results in the paper *Mixture of experts Wishart model for covariance data: applications to Stock markets and Cancer drug screening* by 'Mai T. T., Zhao Z. (2026+)' [arXiv:](https://arxiv.org/abs/).

# Data
## Abstract

The main data contain simulated data and real data. 
The real data is the Cancer Therapeutics Response Portal (CTRP v2, 2015) (Seashore-Ludlow et al., 2015) which is publicly available at 
```diff
https://ocg.cancer.gov/programs/ctd2/data-portal
```

## New R package
We developed an R package `moewishart` which is available on the [GitHub](https://github.com/zhizuio/moewishart/) for implementing our methods. 

## Availability

No restrictions.

# Code 
## Abstract

We are including all of the code that will enable reproducing our simulation results.

# Instructions for Use
## Reproducibility of simulations

Tables 1-2 and Figures 1-2 can be reproduced through the provided code. The general steps are:

1. Install and load our new R-package [moewishart](https://github.com/zhizuio/moewishart/).
2. Load scripts `moewishart_errors.R` and `moewishartX_errors.R`.
3. Run scripts `simN200_p2.R`, `simN500_p2.R` and `simN1000_p2.R`, which will save simulation results in a user-defined folder.
4. Run script `sim_Figure1A.R` and produce Figure 1A. Similarly Figure 1B can be produced by running `simN200_p8.R`, `simN500_p8.R`, `simN1000_p8.R` and slightly modified `sim_Figure1A.R`.
5. Run scripts `simX_N200_p2.R`, `simX_N500_p2.R` and `simX_N1000_p2.R`, which will save simulation results in a user-defined folder.
6. Run script `sim_Figure2A.R` and produce Figure 2A. Similarly Figure 2B can be produced.

## Reproducibility of real data analysis

Table 3 and Figures 4-5 can be reproduced through the provided code. The general steps are:

1. Install and load our new R-package [moewishart](https://github.com/zhizuio/moewishart/).
2. Download the preprocessed CTRP data `CTRP_viability_S_list.RData` and `metaData_selected.csv` from our [GitHub](https://github.com/zhizuio/moewishart/tree/main/simulation_doc), and raw CTRP data `v20.meta.per_compound.csv` and `v20.data.curves_post_qc.txt` from their official [data portal](https://ocg.cancer.gov/programs/ctd2/data-portal).
3. Run script `CTRP_EMs.R` to produce results from **EM** and **EM-MoE** methods.
4. Run script `CTRP_MM.R` to produce results from **Bayes** method.
5. Run script `CTRP_MoE.R` to produce results from **Bayes-MoE** method.
