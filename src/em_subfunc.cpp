
#include <string>
#include <algorithm>
#include <cmath>

#include "em_subfunc.h"

using namespace Rcpp;
using namespace arma;

/*
static std::string to_lower_str(const std::string &s_in)
{
    std::string s = s_in;
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c)
    {
        return std::tolower(c);
    });
    return s;
}
*/
// Multivariate log-gamma: log Gamma_p(a) = (p(p-1)/4) log(pi) + sum_{j=1..p} lgamma(a - (j-1)/2)
double log_multivariate_gamma(int p, double a)
{
    double acc = 0.0;
    for (int j = 0; j < p; ++j)
    {
        acc += std::lgamma(a - (double)j/2.0);
    }
    acc += ( (double)p * (p-1) / 4.0 ) * std::log(datum::pi);
    return acc;
}

// multivariate digamma: psi_p(a) = sum_{j=0..p-1} digamma(a - j/2)
double multivariate_digamma(int p, double a)
{
    double s = 0.0;
    for (int j = 0; j < p; ++j) s += R::digamma(a - (double)j/2.0);
    return s;
}

// multivariate trigamma: psi_p'(a) = sum_{j=0..p-1} trigamma(a - j/2)
double multivariate_trigamma(int p, double a)
{
    double s = 0.0;
    for (int j = 0; j < p; ++j) s += R::trigamma(a - (double)j/2.0);
    return s;
}

// Compute log-density of Wishart_p(S | nu, Sigma)
// Using standard parametrization: density ~ |S|^{(nu-p-1)/2} exp(-0.5 tr(Sigma^{-1} S)) / (2^{nu p/2} |Sigma|^{nu/2} Gamma_p(nu/2))
double log_wishart_density(const arma::mat &S, int p, double nu, const arma::mat &Sigma)
{
    double tol_jitter = 1.0e-10;
    // Ensure SPD Sigma for chol
    arma::mat L;
    bool chol_ok = arma::chol(L, Sigma + tol_jitter * arma::eye<arma::mat>(p,p));
    if (!chol_ok)
    {
        double jitter = 1e-6;
        chol_ok = arma::chol(L, Sigma + jitter * arma::eye<arma::mat>(p,p));
        if (!chol_ok)
        {
            Rcpp::stop("Sigma is not positive definite even after jittering.");
        }
    }
    double logdetSig = 2.0 * arma::sum(arma::log(L.diag()));
    arma::mat Sigma_inv = arma::inv_sympd(Sigma + tol_jitter * arma::eye<arma::mat>(p,p));

    // log|S|
    arma::mat LS;
    double logdetS;
    if (arma::chol(LS, S + tol_jitter * arma::eye<arma::mat>(p,p)))
    {
        logdetS = 2.0 * arma::sum(arma::log(LS.diag()));
    }
    else
    {
        arma::vec svals = arma::svd(S);
        double eps = 1e-12;
        logdetS = 0.0;
        for (int j = 0; j < (int)svals.n_elem; ++j) logdetS += std::log(std::max(svals[j], eps));
    }

    double term1 = ((nu - (double)p - 1.0) / 2.0) * logdetS;
    double term2 = - (nu / 2.0) * logdetSig;
    double trace_term = - 0.5 * arma::trace(Sigma_inv * S);
    double term3 = - ( (nu * (double)p) / 2.0 ) * std::log(2.0);
    double term4 = - log_multivariate_gamma(p, nu / 2.0);

    double logdens = term1 + term2 + trace_term + term3 + term4;
    return logdens;
}

