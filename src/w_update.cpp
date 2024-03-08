#include "RcppArmadillo.h"
#include "CWVS.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List w_update(arma::vec y,
                    arma::mat x,
                    arma::mat z,
                    arma::vec off_set,
                    arma::vec tri_als,
                    int likelihood_indicator,
                    int r,
                    arma::vec beta_old,
                    arma::vec gamma_old,
                    double A11_old,
                    arma::vec delta1_old){

int n = y.size();
  
arma::vec mean_w = off_set +
                   x*beta_old + 
                   z*(gamma_old%(A11_old*delta1_old));

arma::vec input0 = tri_als;
arma::vec input2 = (r + y);

arma::vec w(n); w.fill(0.00);
arma::vec gamma_l(n); gamma_l.fill(0.00);

if(likelihood_indicator == 0){
  
  w = rcpp_pgdraw(input0,
                  mean_w);
  gamma_l = (y - 0.50*tri_als)/w;
  
  } 

if(likelihood_indicator == 2){
  
  w = rcpp_pgdraw(input2,
                  mean_w);
  gamma_l = 0.50*(y - r)/w;
  
  }

return Rcpp::List::create(Rcpp::Named("w") = w,
                          Rcpp::Named("gamma_l") = gamma_l);

}
































































