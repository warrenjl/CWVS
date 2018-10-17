#include "RcppArmadillo.h"
#include "CWVS.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List w_update(arma::vec y,
                    arma::mat x,
                    arma::mat z,
                    arma::vec beta_old,
                    arma::vec gamma_old,
                    double A11_old,
                    arma::vec delta1_old){

arma::vec mean_w = x*beta_old + 
                   z*(gamma_old%(A11_old*delta1_old));

arma::vec w = rcpp_pgdraw(1.0,
                          mean_w);

arma::vec gamma_l = (y - 0.50)/w;

return Rcpp::List::create(Rcpp::Named("w") = w,
                          Rcpp::Named("gamma_l") = gamma_l);

}
































































