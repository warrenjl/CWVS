#include "RcppArmadillo.h"
#include "CWVS.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List A21_update(double A21_old,
                      arma::vec gamma_star,
                      arma::vec delta1,
                      double A22,
                      arma::vec delta2,
                      double sigma2_A,
                      double metrop_var_A21,
                      double acctot_A21){

/*Second*/
arma::vec mean_piece_old = gamma_star - A21_old*delta1 - A22*delta2;

double second = -0.50*dot(mean_piece_old, mean_piece_old) - 
                0.50*(1/sigma2_A)*(A21_old*A21_old);

/*First*/
double A21 = R::rnorm(A21_old, 
                      sqrt(metrop_var_A21));
arma::vec mean_piece = gamma_star - A21*delta1 - A22*delta2;

double first = -0.50*dot(mean_piece, mean_piece) - 
               0.50*(1/sigma2_A)*(A21*A21);

/*Decision*/
double ratio = exp(first - second);   
double acc = 1;
if(ratio < R::runif(0, 1)){
  A21 = A21_old;
  acc = 0;
  }
acctot_A21 = acctot_A21 + 
             acc;

return Rcpp::List::create(Rcpp::Named("A21") = A21,
                          Rcpp::Named("acctot_A21") = acctot_A21);

}



