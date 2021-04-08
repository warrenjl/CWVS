#include "RcppArmadillo.h"
#include "CWVS.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double sigma2_epsilon_update(double a_sigma2_epsilon,
                             double b_sigma2_epsilon,
                             arma::vec y,
                             arma::mat x,
                             arma::mat z,
                             arma::vec beta_old,
                             arma::vec gamma_old,
                             double A11_old,
                             arma::vec delta1_old){
  
int n = y.size();
double a_sigma2_epsilon_update = 0.50*n + 
                                 a_sigma2_epsilon;

double b_sigma2_epsilon_update = 0.50*dot((y - x*beta_old - z*(gamma_old%(A11_old*delta1_old))), (y - x*beta_old - z*(gamma_old%(A11_old*delta1_old)))) + 
                                 b_sigma2_epsilon;

double sigma2_epsilon = 1/R::rgamma(a_sigma2_epsilon_update,
                                    (1.00/b_sigma2_epsilon_update));

return(sigma2_epsilon);

}





