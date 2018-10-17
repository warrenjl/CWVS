#include "RcppArmadillo.h"
#include "CWVS.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

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
                        arma::mat corr_inv1){

int p_z = z.n_cols;
int n = w.size();

arma::mat gamma_mat(n, p_z);
for(int j = 0; j < n; ++j){
   gamma_mat.row(j) = trans(gamma);
   }

arma::mat w_mat(n, p_z);
for(int j = 0; j < p_z; ++j){
  w_mat.col(j) = w;
  }

arma::mat z_gamma = A11_old*(gamma_mat%z);

arma::mat z_gamma_trans = trans(z_gamma);

arma::mat cov_delta1 = inv_sympd(z_gamma_trans*(w_mat%z_gamma) + 
                                 A21_old*A21_old*eye(p_z, p_z) + 
                                 corr_inv1);

arma::vec mean_delta1 = cov_delta1*(z_gamma_trans*(w%(gamma_l - x*beta)) + 
                        A21_old*(gamma_star - A22_old*delta2_old));

arma::mat ind_norms = arma::randn(1, p_z);
arma::vec delta1 = mean_delta1 +                   
                   trans(ind_norms*arma::chol(cov_delta1));

return(delta1);

}

