#include "RcppArmadillo.h"
#include "CWVS.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec delta2_update(arma::vec gamma_star,
                        arma::vec delta1,
                        double A21_old,
                        double A22_old,
                        arma::mat corr_inv2){
  
double p_z = gamma_star.size();
  
arma::mat cov_delta2 = inv_sympd(A22_old*A22_old*eye(p_z, p_z) + 
                                 corr_inv2);

arma::vec mean_delta2 = cov_delta2*(A22_old*(gamma_star - A21_old*delta1));

arma::mat ind_norms = arma::randn(1, p_z);
arma::vec delta2 = mean_delta2 + 
                   trans(ind_norms*arma::chol(cov_delta2));

return(delta2);

}






