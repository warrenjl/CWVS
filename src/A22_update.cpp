#include "RcppArmadillo.h"
#include "CWVS.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List A22_update(double A22_old,
                      arma::vec gamma_star,
                      arma::vec delta1,
                      double A21_old,
                      arma::vec delta2,
                      double sigma2_A,
                      double metrop_var_A22_trans,
                      int acctot_A22_trans){

/*Second*/
double A22_trans_old = log(A22_old);
  
arma::vec mean_piece_old = gamma_star - 
                           A21_old*delta1 - 
                           A22_old*delta2;

double second = -0.50*dot(mean_piece_old, mean_piece_old) - 
                0.50*(1.00/sigma2_A)*(A22_trans_old*A22_trans_old);

/*First*/
double A22_trans = R::rnorm(A22_trans_old, 
                            sqrt(metrop_var_A22_trans));
double A22 = exp(A22_trans);
arma::vec mean_piece = gamma_star - 
                       A21_old*delta1 - 
                       A22*delta2;

double first = -0.50*dot(mean_piece, mean_piece) - 
               0.50*(1.00/sigma2_A)*(A22_trans*A22_trans);

/*Decision*/
double ratio = exp(first - second);   
int acc = 1;
if(ratio < R::runif(0.00, 1.00)){
  A22 = A22_old;
  acc = 0;
  }
acctot_A22_trans = acctot_A22_trans + 
                   acc;

return Rcpp::List::create(Rcpp::Named("A22") = A22,
                          Rcpp::Named("acctot_A22_trans") = acctot_A22_trans);

}



