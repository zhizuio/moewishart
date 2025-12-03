#include <RcppArmadillo.h>
#include <string>
#include <algorithm>
#include <cmath>

#include "em_subfunc.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//' EM for mixture of Wishart distributions with automatic nu_k estimation
//'
//' @name moewishart_em_cpp

//' @param S_list  list of n symmetric positive definite p x p matrices
//' @param K       number of components
//' @param nu_in   scalar or vector length K (initial nu). If NULL, defaults to p+2 per component
//' @param max_em  maximum EM iterations
//' @param max_nr  maximum Newton iterations for nu update
//' @param tol_em  relative tolerance on log-likelihood improvement
//' @param verbose print progress, the.tien.mai@fhi.no
//'
// [[Rcpp::export]]
Rcpp::List moewishart_em_cpp(Rcpp::List S_list,
                             int K = 2,
                             Rcpp::Nullable< Rcpp::NumericVector > nu_in = R_NilValue,
                             int max_em = 200,
                             int max_nr = 50,
                             double tol_em = 1e-8,
                             bool verbose = true)
{
    int n = (int)S_list.size();
    if (n == 0) Rcpp::stop("S_list is empty.");
    arma::mat S0 = as<arma::mat>(S_list[0]);
    int p = S0.n_rows;
    if (S0.n_cols != p) Rcpp::stop("S matrices must be square.");
    // Read/initialize nu vector
    arma::vec nu(K);
    if (nu_in.isNull())
    {
        nu.fill((double)p + 2.0);
    }
    else
    {
        Rcpp::NumericVector nv = as<Rcpp::NumericVector>(nu_in);
        if ((int)nv.size() == 1)
        {
            nu.fill((double)nv[0]);
        }
        else if ((int)nv.size() == K)
        {
            for (int k = 0; k < K; ++k) nu[k] = nv[k];
        }
        else
        {
            Rcpp::stop("nu_in must be scalar or length K.");
        }
    }

    // Convert to arma::field for convenience
    arma::field<arma::mat> Sfield(n);
    arma::vec logdetS_vec(n);
    for (int i = 0; i < n; ++i)
    {
        arma::mat Si = as<arma::mat>(S_list[i]);
        if (Si.n_rows != p || Si.n_cols != p) Rcpp::stop("All S matrices must be same dimension p x p.");
        Sfield(i) = Si;
        // compute logdetS once
        arma::mat LS;
        if (arma::chol(LS, Si + 1e-12 * arma::eye<arma::mat>(p,p)))
        {
            logdetS_vec[i] = 2.0 * arma::sum(arma::log(LS.diag()));
        }
        else
        {
            arma::vec svals = arma::svd(Si);
            double eps = 1e-12;
            double acc = 0.0;
            for (int j = 0; j < (int)svals.n_elem; ++j) acc += std::log(std::max(svals[j], eps));
            logdetS_vec[i] = acc;
        }
    }

    // initialize Sigma_k and pi
    std::vector<arma::mat> Sigma(K);
    arma::vec pi = arma::ones<arma::vec>(K) / (double)K;
    arma::mat meanS = arma::zeros<arma::mat>(p,p);
    for (int i = 0; i < n; ++i) meanS += Sfield(i);
    meanS /= (double)n;
    for (int k = 0; k < K; ++k)
    {
        Sigma[k] = meanS + 0.01 * arma::randn<arma::mat>(p,p);
        Sigma[k] = 0.5 * (Sigma[k] + Sigma[k].t()) + 0.05 * arma::eye<arma::mat>(p,p);
        arma::mat Ltest;
        if (!arma::chol(Ltest, Sigma[k]))
        {
            Sigma[k] += 1e-6 * arma::eye<arma::mat>(p,p);
        }
    }

    arma::mat resp = arma::zeros<arma::mat>(n, K);
    double prev_loglik = -datum::inf;
    double loglik = -datum::inf;

    if (verbose) Rcout << "Mixture-Wishart init: n=" << n << " p=" << p << " K=" << K << "\n";

    for (int em = 0; em < max_em; ++em)
    {
        // E-step: compute log weights
        arma::mat logw = arma::zeros<arma::mat>(n, K);
        for (int k = 0; k < K; ++k)
        {
            double nuk = nu[k];
            for (int i = 0; i < n; ++i)
            {
                double ld = log_wishart_density(Sfield(i), p, nuk, Sigma[k]);
                logw(i, k) = std::log(std::max(pi[k], 1e-300)) + ld;
            }
        }
        // normalize with log-sum-exp
        for (int i = 0; i < n; ++i)
        {
            double m = logw.row(i).max();
            arma::rowvec exps = arma::exp(logw.row(i) - m);
            double s = arma::accu(exps);
            if (s == 0.0) resp.row(i).fill(1.0 / (double)K);
            else resp.row(i) = exps / s;
        }

        // compute loglik
        double new_loglik = 0.0;
        for (int i = 0; i < n; ++i)
        {
            double m = logw.row(i).max();
            double s = arma::accu(arma::exp(logw.row(i) - m));
            new_loglik += m + std::log(s);
        }
        loglik = new_loglik;
        if (verbose) Rcout << "EM iter " << em << ", loglik=" << loglik << "\n";

        // convergence check
        if (em > 0)
        {
            double rel = std::abs((loglik - prev_loglik) / std::max(1.0, std::abs(prev_loglik)));
            if (rel < tol_em && rel > 0.)
            {
                if (verbose) Rcout << "Converged at iter " << em << ", loglik = " << loglik << ", rel = " << rel << "\n";
                break;
            }
        }
        prev_loglik = loglik;

        // M-step: update pi and Sigma and nu
        arma::vec nk = arma::sum(resp, 0).t(); // K-length
        for (int k = 0; k < K; ++k)
        {
            double nk_k = nk[k];
            if (nk_k <= 1e-12)
            {
                // tiny weight: leave parameters unchanged (or small regularization)
                if (verbose) Rcout << "Component " << k << " has tiny nk=" << nk_k << " ; skipping update\n";
                continue;
            }
            // update mixing weight
            pi[k] = nk_k / (double)n;

            // accumulate weighted sum of S_i and weighted avg logdetS
            arma::mat acc = arma::zeros<arma::mat>(p,p);
            double acc_logdet = 0.0;
            for (int i = 0; i < n; ++i)
            {
                acc += resp(i,k) * Sfield(i);
                acc_logdet += resp(i,k) * logdetS_vec[i];
            }
            double mean_logdetS = acc_logdet / nk_k;

            // update Sigma_k: Sigma_k = (1 / (nu_k * nk_k)) * sum_i r_ik * S_i
            arma::mat Sigk = acc / (nu[k] * nk_k);
            Sigk = 0.5 * (Sigk + Sigk.t());

            arma::mat Ltest;
            double jitter = 1e-8;
            bool ok = arma::chol(Ltest, Sigk);

            for (int tries = 0; tries < 10 && !ok; ++tries)
            {
                Sigk += jitter * arma::eye<arma::mat>(p,p);
                ok = arma::chol(Ltest, Sigk);
                jitter *= 10.0;
            }

            Sigma[k] = Sigk;

            // compute logdet Sigma_k
            arma::mat L;
            // bool chol_ok = arma::chol(L, Sigma[k] + 1e-12 * arma::eye<arma::mat>(p,p));
            double logdetSigma = 2.0 * arma::sum(arma::log(L.diag()));

            // Now solve for a = nu/2 via Newton-Raphson on:
            // multivariate_digamma(p, a) = mean_logdetS - logdetSigma - p*log 2
            double rhs = mean_logdetS - logdetSigma - (double)p * std::log(2.0);

            // initial a
            double a = std::max((nu[k] / 2.0), (double)( (p + 1) / 2.0 + 1e-6 )); // safe start
            // lower bound: a > (p-1)/2
            double a_min = (double)( (p - 1) ) / 2.0 + 1e-6;

            bool converged = false;
            for (int it = 0; it < max_nr; ++it)
            {
                // compute g(a) = psi_p(a) - rhs
                double g = multivariate_digamma(p, a) - rhs;
                double gprime = multivariate_trigamma(p, a);
                if (gprime <= 0)   // numerically bad
                {
                    break;
                }
                double step = g / gprime;
                // damped update
                double a_new = a - step;
                // enforce lower bound
                if (a_new <= a_min)
                {
                    a_new = 0.5 * (a + a_min);
                }
                // limit step size to avoid wild jumps
                double max_step = std::max(0.5 * a, 1.0);
                if (std::abs(a_new - a) > max_step)
                {
                    if (a_new > a) a_new = a + max_step;
                    else a_new = a - max_step;
                }
                // compute residual for convergence
                double g_new = multivariate_digamma(p, a_new) - rhs;
                if (std::abs(g_new) < 1e-8)
                {
                    a = a_new;
                    converged = true;
                    break;
                }
                // accept update if smaller |g|
                if (std::abs(g_new) < std::abs(g))
                {
                    a = a_new;
                }
                else
                {
                    // reduce step (damping)
                    a = 0.5 * (a + a_new);
                }
                // ensure a remains > a_min
                if (a <= a_min) a = a_min + 1e-6;
            } // end NR

            if (!converged)
            {
                if (verbose) Rcout << "Warning: nu update for component " << k << " did not fully converge. Keeping previous nu.\n";
                // keep previous nu[k]
            }
            else
            {
                double new_nu = 2.0 * a;
                // enforce reasonable bounds, e.g., nu >= p (or >= p+eps)
                double nu_min = (double)p + 1e-3;
                if (new_nu < nu_min) new_nu = nu_min;
                // optional upper bound to avoid overflow
                double nu_max = 1e6;
                if (new_nu > nu_max) new_nu = nu_max;
                if (verbose) Rcout << "Updated nu[" << k << "] from " << nu[k] << " to " << new_nu << "\n";
                nu[k] = new_nu;
            }
        } // end components loop
    } // end EM loop

    // prepare output
    Rcpp::List Sigma_out(K);
    for (int k = 0; k < K; ++k) Sigma_out[k] = Sigma[k];

    return Rcpp::List::create(
               _["Sigma"] = Sigma_out,
               _["pi"] = pi,
               _["nu"] = nu,
               _["loglik"] = loglik,
               _["responsibilities"] = resp
           );
}
