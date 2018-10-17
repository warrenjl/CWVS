#include "RcppArmadillo.h"
#include "CWVS.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List phi_update(double phi_old,
                      arma::vec delta,
                      Rcpp::List temporal_corr_info,
                      double alpha_phi,
                      double beta_phi,
                      double metrop_var_phi_trans,
                      double acctot_phi_trans){

/*Second*/
Rcpp::List temporal_corr_info_old = temporal_corr_info;
arma::mat corr_inv_old = temporal_corr_info_old[0];
double log_deter_old = temporal_corr_info_old[1];
double phi_trans_old = log(phi_old);

double second = -0.50*log_deter_old - 
                0.50*dot(delta, (corr_inv_old*delta)) + 
                alpha_phi*phi_trans_old -
                beta_phi*exp(phi_trans_old);

/*First*/
double phi_trans = R::rnorm(phi_trans_old, 
                            sqrt(metrop_var_phi_trans));
double phi = exp(phi_trans);
temporal_corr_info = temporal_corr_fun(delta.size(), phi);
arma::mat corr_inv = temporal_corr_info[0];
double log_deter = temporal_corr_info[1];

double first = -0.50*log_deter - 
               0.50*dot(delta, (corr_inv*delta)) + 
               alpha_phi*phi_trans -
               beta_phi*exp(phi_trans);

/*Decision*/
double ratio = exp(first - second);   
double acc = 1;
if(ratio < R::runif(0, 1)){
  phi = phi_old;
  temporal_corr_info = temporal_corr_info_old;
  acc = 0;
  }
acctot_phi_trans = acctot_phi_trans + 
                   acc;

return Rcpp::List::create(Rcpp::Named("phi") = phi,
                          Rcpp::Named("acctot_phi_trans") = acctot_phi_trans,
                          Rcpp::Named("temporal_corr_info") = temporal_corr_info);

}
                 
  
