#include "RcppArmadillo.h"
#include "CWVS.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

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
                       arma::vec delta2_old){
 
arma::vec pieces(2); pieces.fill(0);
arma::vec log_pi(2); log_pi.fill(0);
arma::vec gamma = gamma_old;
arma::vec alpha = gamma%(A11_old*delta1_old);
arma::vec eta = A21_old*delta1_old +
                A22_old*delta2_old;
arma::vec probs(2);
int p_z = z.n_cols;

for(int j = 0; j < p_z; ++j){

   pieces.fill(0);
   log_pi(0) = log(1 - R::pnorm(eta(j),
                                0,
                                1,
                                true,
                                false));
   log_pi(1) = log(R::pnorm(eta(j),
                            0,
                            1,
                            true,
                            false));
   
   for(int k = 0; k < 2; ++k){
      gamma(j) = k;
      alpha(j) = gamma(j)*(A11_old*delta1_old(j));
      pieces(k) = -0.50*dot((gamma_l - x*beta - z*alpha), w%(gamma_l - x*beta - z*alpha)) +
                  log_pi(k);
      }

   probs.fill(0);

   for(int k = 0; k < 2; ++k){
      probs(k) = 1/(sum(exp(pieces - pieces(k))));
  
      if(arma::is_finite(probs(k)) == 0){
        probs(k) = 0;  /*Computational Correction*/
        }
      }

   gamma(j) = as<double>(Rcpp::rbinom(1,
                                      1,
                                      probs(1)));
   alpha(j) = gamma(j)*(A11_old*delta1_old(j));
  
   }

return gamma;

}



