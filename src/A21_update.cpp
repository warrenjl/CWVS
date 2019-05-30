#include "RcppArmadillo.h"
#include "CWVS.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double A21_update(arma::vec gamma_star,
                  arma::vec delta1,
                  double A22,
                  arma::vec delta2,
                  double sigma2_A){
  
double A21_var = 1/(sum(delta1%delta1) + (1.00/sigma2_A));
double A21_mean = A21_var*sum(delta1%(gamma_star - A22*delta2));

double A21 = R::rnorm(A21_mean, 
                      sqrt(A21_var));

return A21;

}



