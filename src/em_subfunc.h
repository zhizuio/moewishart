#ifndef EM_SUBFUNC_H
#define EM_SUBFUNC_H

#include <RcppArmadillo.h>
#include <string>


static std::string to_lower_str(const std::string &s_in);

// Multivariate log-gamma: log Gamma_p(a) = (p(p-1)/4) log(pi) + sum_{j=1..p} lgamma(a - (j-1)/2)
double log_multivariate_gamma(int p, double a);

// multivariate digamma: psi_p(a) = sum_{j=0..p-1} digamma(a - j/2)
double multivariate_digamma(int p, double a);

// multivariate trigamma: psi_p'(a) = sum_{j=0..p-1} trigamma(a - j/2)
double multivariate_trigamma(int p, double a);

// Compute log-density of Wishart_p(S | nu, Sigma)
// Using standard parametrization: density ~ |S|^{(nu-p-1)/2} exp(-0.5 tr(Sigma^{-1} S)) / (2^{nu p/2} |Sigma|^{nu/2} Gamma_p(nu/2))
double log_wishart_density(const arma::mat &S, int p, double nu, const arma::mat &Sigma);

#endif