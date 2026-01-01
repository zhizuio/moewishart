// Main function for full Bayesian


#include <cmath>
#include <vector>
// #include <ctime>

#include "bayes_subfunc.h"

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Main function implemented in C++
//'
//' @name moewishart_cpp TBA
//'
//' @param S_list TBA
//' @param K TBA
//' @param niter TBA
//' @param burnin TBA
//' @param thin TBA
//' @param alpha_in TBA
//' @param nu0_in TBA
//' @param Psi0_in TBA
//' @param init_nu_in TBA
//' @param marginal_z TBA
//' @param sample_nu TBA
//' @param nu_prior_a TBA
//' @param nu_prior_b TBA
//' @param mh_sigma TBA
//' @param verbose TBA
//'
// [[Rcpp::export]]
Rcpp::List moewishart_cpp(Rcpp::List S_list,
                          int K,
                          int niter = 3000,
                          int burnin = 1000,
                          int thin = 1,
                          Rcpp::Nullable<Rcpp::NumericVector> alpha_in = R_NilValue,
                          Rcpp::Nullable<double> nu0_in = R_NilValue,
                          Rcpp::Nullable<Rcpp::NumericMatrix> Psi0_in = R_NilValue,
                          Rcpp::Nullable<Rcpp::NumericVector> init_nu_in = R_NilValue,
                          bool marginal_z = true,
                          bool sample_nu = true,
                          double nu_prior_a = 2.0,
                          double nu_prior_b = 0.1,
                          double mh_sigma = 0.1,
                          bool verbose = true)
{

    int n = S_list.size();
    // determine p from first matrix
    Rcpp::NumericMatrix first = Rcpp::as<Rcpp::NumericMatrix>(S_list[0]);
    int p = first.nrow();

    // alpha
    arma::vec alpha(K);
    if (alpha_in.isNotNull())
    {
        Rcpp::NumericVector av = alpha_in.get();
        for (int k = 0; k < K; ++k) alpha[k] = av[k % av.size()];
    }
    else
    {
        alpha.fill(1.0);
    }

    // CHANGED: nu0 to double
    int nu0;
    if (nu0_in.isNotNull())
    {
        nu0 = Rcpp::as<int>(nu0_in.get());
    }
    else
    {
        nu0 = p + 2;
    }

    arma::mat Psi0(p, p, arma::fill::eye);
    if (Psi0_in.isNotNull())
    {
        Rcpp::NumericMatrix P = Psi0_in.get();
        for (int i = 0; i < p; ++i)
            for (int j = 0; j < p; ++j)
                Psi0(i, j) = P(i, j);
    }

    arma::vec nu_k(K);
    if (init_nu_in.isNotNull())
    {
        Rcpp::NumericVector nv = init_nu_in.get();
        for (int k = 0; k < K; ++k) nu_k[k] = nv[k % nv.size()];
    }
    else
    {
        for (int k = 0; k < K; ++k) nu_k[k] = (double)(p + 2);
    }

    // Vectorize S_list into an n x p^2 arma::mat
    arma::mat S_mat(n, p * p);
    arma::vec log_det_S(n);
    for (int i = 0; i < n; ++i)
    {
        Rcpp::NumericMatrix Si = Rcpp::as<Rcpp::NumericMatrix>(S_list[i]);
        arma::mat Sia(p, p);
        for (int r = 0; r < p; ++r)
            for (int c = 0; c < p; ++c)
                Sia(r, c) = Si(r, c);

        // store as row-major vectorized (column-major of arma: use as_col)
        arma::vec v = arma::vectorise(Sia);
        S_mat.row(i) = v.t();

        double sign;
        arma::log_det(log_det_S[i], sign, Sia);
    }

    // initialize z randomly (uniform)
    arma::ivec z(n);
    for (int i = 0; i < n; ++i) z[i] = (int) (R::unif_rand() * K) + 1; // 1..K

    arma::vec pi_k = rdirichlet_cpp(alpha);

    // Sigma_k: 3D array implemented as vector of arma::mat
    std::vector<arma::mat> Sigma_k(K);
    for (int k = 0; k < K; ++k) Sigma_k[k] = Psi0 / (nu0 - p - 1.0);

    // quick initial update of Sigma_k using simple cluster sums
    for (int k = 0; k < K; ++k)
    {
        arma::uvec idx = arma::find(z == (k + 1));
        if (idx.n_elem > 0)
        {
            arma::rowvec S_sum_vec = arma::sum(S_mat.rows(idx));
            arma::mat S_sum = arma::mat(S_sum_vec.t());
            S_sum.reshape(p, p);
            double nk = (double) idx.n_elem;
            // Note: nu_k is usually parameter of S, nu0 is prior for Sigma.
            // Approximate start:
            Sigma_k[k] = (Psi0 + S_sum) / (nu0 + nk * nu_k[k] - p - 1.0);
        }
    }

    int nsave = std::max(1, (niter - burnin) / thin);
    arma::mat out_pi(nsave, K);
    arma::mat out_nu(nsave, K);
    Rcpp::List out_Sigma(nsave);
    arma::imat out_z(nsave, n);
    arma::vec logliks(niter);
    int iter_save = 0;

    arma::mat logpost(n, K);


    arma::vec counts(K, arma::fill::zeros);
    for (int i = 0; i < n; ++i) counts[z[i] - 1] += 1.0;

    // CHANGED: Use std::clock for timing
    // std::clock_t start_time = std::clock();

    for (int iter = 0; iter < niter; ++iter)
    {
        // --- Step 1: log posterior for each cluster ---
        for (int k = 0; k < K; ++k)
        {
            arma::mat Sig = Sigma_k[k];
            arma::mat cholSig;
            bool ok = arma::chol(cholSig, Sig);
            if (!ok)
            {
                arma::mat Sig2 = Sig + arma::eye<arma::mat>(p, p) * 1e-6;
                arma::chol(cholSig, Sig2);
            }
            double log_det_Sig = 2.0 * arma::sum(arma::log(cholSig.diag()));
            arma::mat Sig_inv = arma::inv(arma::trimatu(cholSig)).t() * arma::inv(arma::trimatu(cholSig));

            // tr(Sig_inv * S_i) for all i: S_mat * vec(Sig_inv)
            arma::vec tr_val = S_mat * arma::vectorise(Sig_inv);

            double nu = nu_k[k];
            // term1: (nu - p - 1)/2 * log|S|
            arma::vec term1 = ((nu - p - 1.0) / 2.0) * log_det_S;
            arma::vec term2 = -0.5 * tr_val;
            double term3_const = - (nu * p / 2.0) * std::log(2.0) - (nu / 2.0) * log_det_Sig;
            double term4_const = - log_multivariate_gamma(nu / 2.0, p);

            for (int i = 0; i < n; ++i)
            {
                if (marginal_z)
                {
                    // using (partially) marginal z by integrating out pi_k
                    logpost(i, k) = std::log(alpha[k] + counts[k] + 1e-300) + term1[i] + term2[i] + term3_const + term4_const;
                }
                else
                {
                    logpost(i, k) = std::log(pi_k[k] + 1e-300) + term1[i] + term2[i] + term3_const + term4_const;
                }
            }
        }

        // Sample z for each i
        for (int i = 0; i < n; ++i)
        {
            arma::rowvec lp = logpost.row(i);
            double m = lp.max();
            arma::rowvec shifted = lp - m;
            arma::rowvec pvec = arma::exp(shifted);
            double s = arma::accu(pvec);
            if (s <= 0) pvec.fill(1.0 / K);
            else pvec /= s;

            // sample from categorical
            double u = R::unif_rand();
            double csum = 0.0;
            int choose = 0;
            for (int k = 0; k < K; ++k)
            {
                csum += pvec[k];
                if (u <= csum)
                {
                    choose = k;
                    break;
                }
            }
            z[i] = choose + 1;
        }

        // Step 2: Update pi
        arma::vec counts(K, arma::fill::zeros);
        for (int i = 0; i < n; ++i) counts[z[i] - 1] += 1.0;
        arma::vec alpha_post = alpha + counts;
        pi_k = rdirichlet_cpp(alpha_post);

        // Step 3: Update Sigma_k
        for (int k = 0; k < K; ++k)
        {
            arma::uvec idx = arma::find(z == (k + 1));
            int nk = idx.n_elem;
            if (nk == 0)
            {
                Sigma_k[k] = sample_invwishart(nu0, Psi0);
            }
            else
            {
                arma::rowvec S_sum_vec = arma::sum(S_mat.rows(idx));
                arma::mat S_sum = arma::mat(S_sum_vec.t());
                S_sum.reshape(p, p);
                // CHANGED: Use double arithmetic for nu_post to match MH sampler
                double nu_post = nu0 + (double)nk * nu_k[k];
                arma::mat Psi_post = Psi0 + S_sum;
                // sample inverse-Wishart with (nu_post, Psi_post)
                Sigma_k[k] = sample_invwishart(nu_post, Psi_post);
            }
        }

        // Step 4: Update nu via random-walk on log-scale (MH)
        if (sample_nu)
        {
            for (int k = 0; k < K; ++k)
            {
                double curr_nu = nu_k[k];
                double prop_log = R::rnorm(std::log(curr_nu), mh_sigma);
                double prop_nu = std::exp(prop_log);

                if (prop_nu > (p - 1.0 + 1e-6))
                {
                    arma::uvec idx = arma::find(z == (k + 1));
                    double ll_old = 0.0, ll_new = 0.0;
                    if (idx.n_elem > 0)
                    {
                        double sum_log_det_S_k = arma::sum(log_det_S.elem(idx));

                        // re-invert current Sigma_k
                        arma::mat Sig = Sigma_k[k];
                        arma::mat cholSig;
                        bool ok = arma::chol(cholSig, Sig);
                        if (!ok) arma::chol(cholSig, Sig + arma::eye<arma::mat>(p, p) * 1e-6);
                        double log_det_Sig = 2.0 * arma::sum(arma::log(cholSig.diag()));
                        arma::mat Sig_inv = arma::inv(arma::trimatu(cholSig)).t() * arma::inv(arma::trimatu(cholSig));

                        arma::rowvec S_sum_vec = arma::sum(S_mat.rows(idx));
                        arma::vec Ssumvec = S_sum_vec.t();
                        double sum_tr_val = arma::dot(arma::vectorise(Sig_inv), Ssumvec);

                        // Lambda to calc LL for a specific nu
                        auto calc_ll_nu = [&](double val_nu)
                        {
                            double term1 = (val_nu - p - 1.0) / 2.0 * sum_log_det_S_k;
                            double term2 = -0.5 * sum_tr_val;
                            double term3 = -((double)idx.n_elem) * ((val_nu * p / 2.0) * std::log(2.0) + (val_nu / 2.0) * log_det_Sig);
                            double term4 = -((double)idx.n_elem) * log_multivariate_gamma(val_nu / 2.0, p);
                            return term1 + term2 + term3 + term4;
                        };

                        ll_old = calc_ll_nu(curr_nu);
                        ll_new = calc_ll_nu(prop_nu);
                    }

                    double lp_old = (nu_prior_a - 1.0) * std::log(curr_nu) - nu_prior_b * curr_nu ; // Jacobian
                    double lp_new = (nu_prior_a - 1.0) * std::log(prop_nu) - nu_prior_b * prop_nu ; // Jacobian

                    double accept_logprob = (ll_new + lp_new) - (ll_old + lp_old);
                    accept_logprob += std::log(prop_nu) - std::log(curr_nu); // add proposal ratio
                    if (std::log(R::runif(0.0, 1.0)) < accept_logprob) nu_k[k] = prop_nu;
                }
            }
        }

        // Log-likelihood approximation
        arma::vec max_l = arma::max(logpost, 1);
        arma::mat row_exps = arma::exp(logpost.each_col() - max_l);
        arma::vec row_sums = arma::sum(row_exps, 1);
        logliks[iter] = arma::accu(max_l + arma::log(row_sums));

        // Save
        if (iter >= burnin && ((iter - burnin) % thin == 0))
        {
            if (iter_save < nsave)
            {
                out_pi.row(iter_save) = pi_k.t();
                out_nu.row(iter_save) = nu_k.t();

                Rcpp::List saveK(K);
                for (int k = 0; k < K; ++k)
                {
                    Rcpp::NumericMatrix M(p, p);
                    // Explicit copy to R matrix to ensure safe storage in List
                    arma::mat Sk = Sigma_k[k];
                    for (int i = 0; i < p; ++i)
                        for (int j = 0; j < p; ++j)
                            M(i, j) = Sk(i, j);
                    saveK[k] = M;
                }
                out_Sigma[iter_save] = saveK;

                for (int i = 0; i < n; ++i) out_z(iter_save, i) = z[i];
                iter_save++;
            }
        }

        if (verbose && (((iter+1) % 500 == 0) || (iter+1 == 0)))
        {
            // CHANGED: Correct timing calculation
            // double elapsed = (double)(std::clock() - start_time) / CLOCKS_PER_SEC;
            // double rate = (elapsed > 0) ? (iter + 1) / elapsed : 0.0;
            Rcpp::Rcout << "Iter " << iter+1 << " | LL=" << logliks[iter] << " | pi=" << pi_k.t();// << "\n";
        }
    }

    // posterior mean of Sigma (average over saved lists)
    std::vector<arma::mat> Sigma_mean(K, arma::zeros<arma::mat>(p, p));
    for (int s = 0; s < std::min(nsave, iter_save); ++s)
    {
        // CHANGED: Safe casting from generic List element
        Rcpp::List saveK = Rcpp::as<Rcpp::List>(out_Sigma[s]);
        for (int k = 0; k < K; ++k)
        {
            Rcpp::NumericMatrix M = Rcpp::as<Rcpp::NumericMatrix>(saveK[k]);
            arma::mat Ma(p, p);
            for (int i = 0; i < p; ++i)
            {
                for (int j = 0; j < p; ++j)
                    Ma(i, j) = M(i, j);
            }
            Sigma_mean[k] += Ma;
        }
    }

    for (int k = 0; k < K; ++k) Sigma_mean[k] /= (double) std::max(1, iter_save);

    // convert Sigma_mean into a single List of K matrices
    Rcpp::List Sigma_mean_list(K);
    for (int k = 0; k < K; ++k)
    {
        Rcpp::NumericMatrix M(p, p);
        for (int i = 0; i < p; ++i)
            for (int j = 0; j < p; ++j)
                M(i, j) = Sigma_mean[k](i, j);
        Sigma_mean_list[k] = M;
    }

    return Rcpp::List::create(
               Rcpp::Named("pi") = out_pi,
               Rcpp::Named("nu") = out_nu,
               Rcpp::Named("Sigma") = out_Sigma,
               Rcpp::Named("z") = out_z,
               Rcpp::Named("sigma_posterior_mean") = Sigma_mean_list,
               Rcpp::Named("loglik") = logliks
           );
}