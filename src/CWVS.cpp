#include "RcppArmadillo.h"
#include "CWVS.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List CWVS(int mcmc_samples,
                arma::vec y,
                arma::mat x,
                arma::mat z,
                int likelihood_indicator,
                double metrop_var_phi1_trans,
                double metrop_var_phi2_trans,
                double metrop_var_A11_trans,
                double metrop_var_A22_trans,
                Rcpp::Nullable<Rcpp::NumericVector> offset = R_NilValue,
                Rcpp::Nullable<Rcpp::NumericVector> trials = R_NilValue,
                Rcpp::Nullable<double> a_r_prior = R_NilValue,
                Rcpp::Nullable<double> b_r_prior = R_NilValue,
                Rcpp::Nullable<double> a_sigma2_epsilon_prior = R_NilValue,
                Rcpp::Nullable<double> b_sigma2_epsilon_prior = R_NilValue,
                Rcpp::Nullable<double> sigma2_beta_prior = R_NilValue,
                Rcpp::Nullable<double> alpha_phi1_prior = R_NilValue,
                Rcpp::Nullable<double> beta_phi1_prior = R_NilValue,
                Rcpp::Nullable<double> alpha_phi2_prior = R_NilValue,
                Rcpp::Nullable<double> beta_phi2_prior = R_NilValue,
                Rcpp::Nullable<double> sigma2_A_prior = R_NilValue,
                Rcpp::Nullable<double> r_init = R_NilValue,
                Rcpp::Nullable<double> sigma2_epsilon_init = R_NilValue,
                Rcpp::Nullable<Rcpp::NumericVector> beta_init = R_NilValue,
                Rcpp::Nullable<Rcpp::NumericVector> gamma_init = R_NilValue,
                Rcpp::Nullable<Rcpp::NumericVector> delta1_init = R_NilValue,
                Rcpp::Nullable<Rcpp::NumericVector> delta2_init = R_NilValue,
                Rcpp::Nullable<double> phi1_init = R_NilValue,
                Rcpp::Nullable<double> phi2_init = R_NilValue,
                Rcpp::Nullable<double> A11_init = R_NilValue,
                Rcpp::Nullable<double> A22_init = R_NilValue,
                Rcpp::Nullable<double> A21_init = R_NilValue){

//Defining Parameters and Quantities of Interest
int n = y.size();
int p_x = x.n_cols;
int p_z = z.n_cols;
double max_time = (p_z - 1);
arma::vec r(mcmc_samples); r.fill(0.00);
arma::vec sigma2_epsilon(mcmc_samples); sigma2_epsilon.fill(0.00);
arma::mat beta(p_x, mcmc_samples); beta.fill(0.00);
arma::mat gamma(p_z, mcmc_samples); gamma.fill(0.00);
arma::mat delta1(p_z, mcmc_samples); delta1.fill(0.00);
arma::mat delta2(p_z, mcmc_samples); delta2.fill(0.00);
arma::vec phi1(mcmc_samples); phi1.fill(0.00);
arma::vec phi2(mcmc_samples); phi2.fill(0.00);
arma::vec A11(mcmc_samples); A11.fill(0.00);
arma::vec A22(mcmc_samples); A22.fill(0.00);
arma::vec A21(mcmc_samples); A21.fill(0.00);
arma::mat alpha(p_z, mcmc_samples); alpha.fill(0.00);
arma::vec neg_two_loglike(mcmc_samples); neg_two_loglike.fill(0.00);

arma::vec off_set(n); off_set.fill(0.00);
if(offset.isNotNull()){
  off_set = Rcpp::as<arma::vec>(offset);
  }

arma::vec tri_als(n); tri_als.fill(1);
if(trials.isNotNull()){
  tri_als = Rcpp::as<arma::vec>(trials);
  }

//Prior Information
int a_r = 1;
if(a_r_prior.isNotNull()){
  a_r = Rcpp::as<int>(a_r_prior);
  }

int b_r = 100;
if(b_r_prior.isNotNull()){
  b_r = Rcpp::as<int>(b_r_prior);
  }

double a_sigma2_epsilon = 0.01;
if(a_sigma2_epsilon_prior.isNotNull()){
  a_sigma2_epsilon = Rcpp::as<double>(a_sigma2_epsilon_prior);
  }

double b_sigma2_epsilon = 0.01;
if(b_sigma2_epsilon_prior.isNotNull()){
  b_sigma2_epsilon = Rcpp::as<double>(b_sigma2_epsilon_prior);
  }

double sigma2_beta = 10000.00;
if(sigma2_beta_prior.isNotNull()){
  sigma2_beta = Rcpp::as<double>(sigma2_beta_prior);
  }

double alpha_phi1 = 1.00;
if(alpha_phi1_prior.isNotNull()){
  alpha_phi1 = Rcpp::as<double>(alpha_phi1_prior);
  }
  
double beta_phi1 = 1.00;
if(beta_phi1_prior.isNotNull()){
  beta_phi1 = Rcpp::as<double>(beta_phi1_prior);
  }

double alpha_phi2 = 1.00;
if(alpha_phi2_prior.isNotNull()){
  alpha_phi2 = Rcpp::as<double>(alpha_phi2_prior);
  }

double beta_phi2 = 1.00;
if(beta_phi2_prior.isNotNull()){
  beta_phi2 = Rcpp::as<double>(beta_phi2_prior);
  }

double sigma2_A = 1.00;
if(sigma2_A_prior.isNotNull()){
  sigma2_A = Rcpp::as<double>(sigma2_A_prior);
  }

//Initial Values
r(0) = b_r;
if(r_init.isNotNull()){
  r(0) = Rcpp::as<int>(r_init);
  }

sigma2_epsilon(0) = 1.00;
if(sigma2_epsilon_init.isNotNull()){
  sigma2_epsilon(0) = Rcpp::as<double>(sigma2_epsilon_init);
  }

beta.col(0).fill(0.00);
if(beta_init.isNotNull()){
  beta.col(0) = Rcpp::as<arma::vec>(beta_init);
  }

gamma.col(0).fill(1.00);
if(gamma_init.isNotNull()){
  gamma.col(0) = Rcpp::as<arma::vec>(gamma_init);
  }

delta1.col(0).fill(0.00);
if(delta1_init.isNotNull()){
  delta1.col(0) = Rcpp::as<arma::vec>(delta1_init);
  }

delta2.col(0).fill(0.00);
if(delta2_init.isNotNull()){
  delta2.col(0) = Rcpp::as<arma::vec>(delta2_init);
  }

phi1(0) = -log(0.05)/max_time;  //Effective range equal to largest temporal distance in dataset (strong temporal correlation)
if(phi1_init.isNotNull()){
  phi1(0) = Rcpp::as<double>(phi1_init);
  }

phi2(0) = -log(0.05)/max_time;  //Effective range equal to largest temporal distance in dataset (strong temporal correlation)
if(phi2_init.isNotNull()){
  phi2(0) = Rcpp::as<double>(phi2_init);
  }

A11(0) = 1.00;
if(A11_init.isNotNull()){
  A11(0) = Rcpp::as<double>(A11_init);
  }

A22(0) = 1.00;
if(A22_init.isNotNull()){
  A22(0) = Rcpp::as<double>(A22_init);
  }

A21(0) = 0.00;
if(A21_init.isNotNull()){
  A21(0) = Rcpp::as<double>(A21_init);
  }

Rcpp::List temporal_corr_info1 = temporal_corr_fun(p_z, phi1(0));
Rcpp::List temporal_corr_info2 = temporal_corr_fun(p_z, phi2(0));
neg_two_loglike(0) = neg_two_loglike_update(y,
                                            x,
                                            z, 
                                            off_set,
                                            likelihood_indicator,
                                            r(0),
                                            sigma2_epsilon(0),
                                            beta.col(0),
                                            gamma.col(0),
                                            A11(0),
                                            delta1.col(0));

//Metropolis Settings
int acctot_phi1_trans = 0;
int acctot_phi2_trans = 0;
int acctot_A11_trans = 0;
int acctot_A22_trans = 0;

//Main Sampling Loop
arma::vec w(y.size()); w.fill(0.00);
arma::vec gamma_l = y;
if(likelihood_indicator == 2){
  
  Rcpp::List w_output = w_update(y,
                                 x,
                                 z,
                                 off_set,
                                 tri_als,
                                 likelihood_indicator,
                                 r(0),
                                 beta.col(0),
                                 gamma.col(0),
                                 A11(0),
                                 delta1.col(0));
  w = Rcpp::as<arma::vec>(w_output[0]);
  gamma_l = Rcpp::as<arma::vec>(w_output[1]);
  
  }

for(int j = 1; j < mcmc_samples; ++j){
  
   if(likelihood_indicator == 1){
    
     //sigma2_epsilon Update
     sigma2_epsilon(j) = sigma2_epsilon_update(a_sigma2_epsilon,
                                               b_sigma2_epsilon,
                                               y,
                                               x,
                                               z,
                                               beta.col(j-1),
                                               gamma.col(j-1),
                                               A11(j-1),
                                               delta1.col(j-1));
     w.fill(1.00/sigma2_epsilon(j));
    
     }
  
   if(likelihood_indicator == 0){
    
     //w Update
     Rcpp::List w_output = w_update(y,
                                    x,
                                    z,
                                    off_set,
                                    tri_als,
                                    likelihood_indicator,
                                    r(j-1),
                                    beta.col(j-1),
                                    gamma.col(j-1),
                                    A11(j-1),
                                    delta1.col(j-1));
  
     w = Rcpp::as<arma::vec>(w_output[0]);
     gamma_l = Rcpp::as<arma::vec>(w_output[1]);
     
     }
  
  //beta Update
  beta.col(j) = beta_update(x, 
                            z,
                            off_set,
                            sigma2_beta,
                            w,
                            gamma_l,
                            gamma.col(j-1),
                            A11(j-1),
                            delta1.col(j-1));
  
  //gamma Update
  gamma.col(j) = gamma_update(x, 
                              z,
                              off_set,
                              w,
                              gamma_l,
                              beta.col(j),
                              gamma.col(j-1),
                              A11(j-1),
                              delta1.col(j-1),
                              A21(j-1),
                              A22(j-1),
                              delta2.col(j-1));
  
  //gamma_star Update
  arma::vec gamma_star = gamma_star_update(gamma.col(j),
                                           delta1.col(j-1),
                                           A21(j-1),
                                           A22(j-1),
                                           delta2.col(j-1));
  
  //delta1 Update
  delta1.col(j) = delta1_update(x,
                                z,
                                off_set,
                                w,
                                gamma_l,
                                beta.col(j),
                                gamma.col(j),
                                gamma_star,
                                A11(j-1),
                                A21(j-1),
                                A22(j-1),
                                delta2.col(j-1),
                                temporal_corr_info1(0));
  
  //phi1 Update
  Rcpp::List phi1_output = phi_update(phi1(j-1),
                                      delta1.col(j),
                                      temporal_corr_info1,
                                      alpha_phi1,
                                      beta_phi1,
                                      metrop_var_phi1_trans,
                                      acctot_phi1_trans);
  
  phi1(j) = Rcpp::as<double>(phi1_output[0]);
  acctot_phi1_trans = phi1_output[1];
  temporal_corr_info1 = phi1_output[2];

  //delta2 Update
  delta2.col(j) = delta2_update(gamma_star,
                                delta1.col(j),
                                A21(j-1),
                                A22(j-1),
                                temporal_corr_info2(0));
  
  //phi2 Update
  Rcpp::List phi2_output = phi_update(phi2(j-1),
                                      delta2.col(j),
                                      temporal_corr_info2,
                                      alpha_phi2,
                                      beta_phi2,
                                      metrop_var_phi2_trans,
                                      acctot_phi2_trans);
  
  phi2(j) = Rcpp::as<double>(phi2_output[0]);
  acctot_phi2_trans = phi2_output[1];
  temporal_corr_info2 = phi2_output[2];
  
  //A11 Update
  Rcpp::List A11_output = A11_update(A11(j-1),
                                     x,
                                     z,
                                     off_set,
                                     w,
                                     gamma_l,
                                     beta.col(j),
                                     gamma.col(j),
                                     delta1.col(j),
                                     sigma2_A,
                                     metrop_var_A11_trans,
                                     acctot_A11_trans);
  
  A11(j) = Rcpp::as<double>(A11_output[0]);
  acctot_A11_trans = A11_output[1];
  
  //A22 Update
  Rcpp::List A22_output = A22_update(A22(j-1),
                                     gamma_star,
                                     delta1.col(j),
                                     A21(j-1),
                                     delta2.col(j),
                                     sigma2_A,
                                     metrop_var_A22_trans,
                                     acctot_A22_trans);
  
  A22(j) = Rcpp::as<double>(A22_output[0]);
  acctot_A22_trans = A22_output[1];
  
  //A21 Update
  A21(j) = A21_update(gamma_star,
                      delta1.col(j),
                      A22(j),
                      delta2.col(j),
                      sigma2_A);
  
  //alpha Update
  alpha.col(j) = gamma.col(j)%(A11(j)*delta1.col(j));
  
  if(likelihood_indicator == 2){
    
    //r Update
    r(j) = r_update(y,
                    x,
                    z,
                    off_set,
                    a_r,
                    b_r,
                    beta.col(j),
                    gamma.col(j),
                    A11(j),
                    delta1.col(j));
    
    //w Update
    Rcpp::List w_output = w_update(y,
                                   x,
                                   z,
                                   off_set,
                                   tri_als,
                                   likelihood_indicator,
                                   r(j),
                                   beta.col(j),
                                   gamma.col(j),
                                   A11(j),
                                   delta1.col(j));
    w = Rcpp::as<arma::vec>(w_output[0]);
    gamma_l = Rcpp::as<arma::vec>(w_output[1]);
    
    }
  
  //neg_two_loglike Update
  neg_two_loglike(j) = neg_two_loglike_update(y,
                                              x,
                                              z, 
                                              off_set,
                                              likelihood_indicator,
                                              r(j),
                                              sigma2_epsilon(j),
                                              beta.col(j),
                                              gamma.col(j),
                                              A11(j),
                                              delta1.col(j));
  
  //Progress
  if((j + 1) % 10 == 0){ 
    Rcpp::checkUserInterrupt();
    }
  
  if(((j + 1) % int(round(mcmc_samples*0.10)) == 0)){
    double completion = round(100*((j + 1)/(double)mcmc_samples));
    Rcpp::Rcout << "Progress: " << completion << "%" << std::endl;
    double accrate_phi1_trans = round(100*(acctot_phi1_trans/(double)j));
    Rcpp::Rcout << "phi1 Acceptance: " << accrate_phi1_trans << "%" << std::endl;
    double accrate_phi2_trans = round(100*(acctot_phi2_trans/(double)j));
    Rcpp::Rcout << "phi2 Acceptance: " << accrate_phi2_trans << "%" << std::endl;
    double accrate_A11_trans = round(100*(acctot_A11_trans/(double)j));
    Rcpp::Rcout << "A11 Acceptance: " << accrate_A11_trans << "%" << std::endl;
    double accrate_A22_trans = round(100*(acctot_A22_trans/(double)j));
    Rcpp::Rcout << "A22 Acceptance: " << accrate_A22_trans << "%" << std::endl;
    Rcpp::Rcout << "********************" << std::endl;
    }
  
  }
                                  
return Rcpp::List::create(Rcpp::Named("sigma2_epsilon") = sigma2_epsilon,
                          Rcpp::Named("beta") = beta,
                          Rcpp::Named("gamma") = gamma,
                          Rcpp::Named("delta1") = delta1,
                          Rcpp::Named("delta2") = delta2,
                          Rcpp::Named("phi1") = phi1,
                          Rcpp::Named("phi2") = phi2,
                          Rcpp::Named("A11") = A11,
                          Rcpp::Named("A22") = A22,
                          Rcpp::Named("A21") = A21,
                          Rcpp::Named("alpha") = alpha,
                          Rcpp::Named("r") = r,
                          Rcpp::Named("neg_two_loglike") = neg_two_loglike,
                          Rcpp::Named("acctot_phi1_trans") = acctot_phi1_trans,
                          Rcpp::Named("acctot_phi2_trans") = acctot_phi2_trans,
                          Rcpp::Named("acctot_A11_trans") = acctot_A11_trans,
                          Rcpp::Named("acctot_A22_trans") = acctot_A22_trans);

}
