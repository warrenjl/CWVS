#ifndef __CWVS__
#define __CWVS__

arma::vec rcpp_pgdraw(double b, 
                      arma::vec c);

Rcpp::List temporal_corr_fun(int p_z,
                             double phi);

Rcpp::List w_update(arma::vec y,
                    arma::mat x,
                    arma::mat z,
                    arma::vec beta_old,
                    arma::vec gamma_old,
                    double A11_old,
                    arma::vec delta1_old);

arma::vec beta_update(arma::mat x, 
                      arma::mat z,
                      double sigma2_beta,
                      arma::vec w,
                      arma::vec gamma_l,
                      arma::vec gamma_old,
                      double A11_old,
                      arma::vec delta1_old);

arma::vec gamma_update(arma::mat x, 
                       arma::mat z,
                       arma::vec w,
                       arma::vec gamma_l,
                       arma::vec beta,
                       arma::vec gamma_old,
                       double A11_old,
                       arma::vec delta1_old,
                       double A21_old,
                       double A22_old,
                       arma::vec delta2_old);

arma::vec gamma_star_update(arma::vec gamma,
                            arma::vec delta1_old,
                            double A21_old,
                            double A22_old,
                            arma::vec delta2_old);

arma::vec delta1_update(arma::mat x,
                        arma::mat z,
                        arma::vec w,
                        arma::vec gamma_l,
                        arma::vec beta,
                        arma::vec gamma,
                        arma::vec gamma_star,
                        double A11_old,
                        double A21_old,
                        double A22_old,
                        arma::vec delta2_old,
                        arma::mat corr_inv1);

Rcpp::List phi_update(double phi_old,
                      arma::vec delta,
                      Rcpp::List temporal_corr_info,
                      double alpha_phi,
                      double beta_phi,
                      double metrop_var_phi_trans,
                      double acctot_phi_trans);

arma::vec delta2_update(arma::vec gamma_star,
                        arma::vec delta1,
                        double A21_old,
                        double A22_old,
                        arma::mat corr_inv2);

Rcpp::List A11_update(double A11_old,
                      arma::mat x,
                      arma::mat z,
                      arma::vec w,
                      arma::vec gamma_l,
                      arma::vec beta,
                      arma::vec gamma,
                      arma::vec delta1,
                      double sigma2_A,
                      double metrop_var_A11_trans,
                      double acctot_A11_trans);

Rcpp::List A22_update(double A22_old,
                      arma::vec gamma_star,
                      arma::vec delta1,
                      double A21_old,
                      arma::vec delta2,
                      double sigma2_A,
                      double metrop_var_A22_trans,
                      double acctot_A22_trans);

Rcpp::List A21_update(double A21_old,
                      arma::vec gamma_star,
                      arma::vec delta1,
                      double A22,
                      arma::vec delta2,
                      double sigma2_A,
                      double metrop_var_A21,
                      double acctot_A21);

double neg_two_loglike_update(arma::vec y,
                              arma::mat x,
                              arma::mat z, 
                              arma::vec beta,
                              arma::vec gamma,
                              double A11,
                              arma::vec delta1);

Rcpp::List CWVS(int mcmc_samples,
                arma::vec y,
                arma::mat x,
                arma::mat z,
                double metrop_var_phi1_trans,
                double metrop_var_phi2_trans,
                double metrop_var_A11_trans,
                double metrop_var_A22_trans,
                double metrop_var_A21,
                Rcpp::Nullable<double> sigma2_beta_prior,
                Rcpp::Nullable<double> alpha_phi1_prior,
                Rcpp::Nullable<double> beta_phi1_prior,
                Rcpp::Nullable<double> alpha_phi2_prior,
                Rcpp::Nullable<double> beta_phi2_prior,
                Rcpp::Nullable<double> sigma2_A_prior,
                Rcpp::Nullable<Rcpp::NumericVector> beta_init,
                Rcpp::Nullable<Rcpp::NumericVector> gamma_init,
                Rcpp::Nullable<Rcpp::NumericVector> delta1_init,
                Rcpp::Nullable<Rcpp::NumericVector> delta2_init,
                Rcpp::Nullable<double> phi1_init,
                Rcpp::Nullable<double> phi2_init,
                Rcpp::Nullable<double> A11_init,
                Rcpp::Nullable<double> A22_init,
                Rcpp::Nullable<double> A21_init); 

#endif // __CWVS__
