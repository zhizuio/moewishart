// subfunctions for full Bayesian

#include <RcppArmadillo.h>
#include "bayes_subfunc.h"


// Log multivariate gamma: log Gamma_p(a)
double log_multivariate_gamma(double a, int p)
{
    double out = (p * (p - 1) / 4.0) * std::log(M_PI);
    for (int j = 0; j < p; ++j)
    {
        out += R::lgammafn(a + (1.0 - (j + 1)) / 2.0);
    }
    return out;
}

// Sample a Dirichlet vector using gamma draws
arma::vec rdirichlet_cpp(const arma::vec &alpha)
{
    int K = alpha.n_elem;
    arma::vec g(K);
    double s = 0.0;
    for (int k = 0; k < K; ++k)
    {
        g[k] = R::rgamma(alpha[k], 1.0);
        s += g[k];
    }
    if (s == 0.0) return arma::ones<arma::vec>(K) / (double)K;
    return g / s;
}

// Sample Wishart(nu, Sigma) using Bartlett decomposition.
arma::mat sample_wishart(double nu, const arma::mat &Sigma)
{
    int p = Sigma.n_rows;
    arma::mat A = arma::zeros<arma::mat>(p, p);
    for (int i = 0; i < p; ++i)
    {
        // diagonal: sqrt of chi-square
        double df = nu - (double)i;
        if (df <= 0) df = 1e-8; // fallback
        A(i, i) = std::sqrt(R::rchisq(df));
        for (int j = 0; j < i; ++j)
        {
            A(i, j) = R::rnorm(0.0, 1.0);
        }
    }

    // Cholesky of Sigma (lower-triangular)
    arma::mat C;
    bool ok = arma::chol(C, Sigma, "lower");
    if (!ok)
    {
        // jitter and retry
        arma::mat S2 = Sigma + arma::eye<arma::mat>(p, p) * 1e-6;
        arma::chol(C, S2, "lower");
    }
    // B = C * A  (A is lower-triangular)
    arma::mat B = C * A;
    arma::mat W = B * B.t();
    return W;
}

// Sample inverse-Wishart(nu, Psi)
arma::mat sample_invwishart(double nu, const arma::mat &Psi)
{
    // Sample W ~ Wishart(nu, inv(Psi)) then return inv(W)
    arma::mat Psi_inv;
    bool ok = arma::inv(Psi_inv, Psi);
    if (!ok) Psi_inv = arma::inv(Psi + arma::eye<arma::mat>(Psi.n_rows, Psi.n_rows) * 1e-6);

    arma::mat W = sample_wishart(nu, Psi_inv);

    arma::mat W_inv;
    ok = arma::inv(W_inv, W);
    if (!ok) W_inv = arma::inv(W + arma::eye<arma::mat>(Psi.n_rows, Psi.n_rows) * 1e-6);
    return W_inv;
}
