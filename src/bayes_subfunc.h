#ifndef BAYES_SUBFUNC_H
#define BAYES_SUBFUNC_H

#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <ctime>

// Log multivariate gamma: log Gamma_p(a)
double log_multivariate_gamma(double a, int p);

// Sample a Dirichlet vector using gamma draws
arma::vec rdirichlet_cpp(const arma::vec &alpha);

// Sample Wishart(nu, Sigma) using Bartlett decomposition.
arma::mat sample_wishart(double nu, const arma::mat &Sigma);

// Sample inverse-Wishart(nu, Psi)
arma::mat sample_invwishart(double nu, const arma::mat &Psi);

#endif