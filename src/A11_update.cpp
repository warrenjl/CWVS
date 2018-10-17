#include "RcppArmadillo.h"
#include "CWVS.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List A11_update(double A11_old,
                      arma::mat x,
                      arma::mat z,
                      arma::vec w,
                      arma::vec gamma_l,
                      arma::vec beta,
                      arma::vec gamma,
                      arma::vec delta1,
                      double sigma2_A,
                      double metrop_var_A11_trans,
                      double acctot_A11_trans){

/*Second*/
double A11_trans_old = log(A11_old);
  
arma::vec mean_piece_old = gamma_l - x*beta - z*(gamma%(A11_old*delta1));

double second = -0.50*dot(mean_piece_old, w%mean_piece_old) - 
                0.50*(1/sigma2_A)*(A11_trans_old*A11_trans_old);

/*First*/
double A11_trans = R::rnorm(A11_trans_old, 
                            sqrt(metrop_var_A11_trans));
double A11 = exp(A11_trans);
arma::vec mean_piece = gamma_l - x*beta - z*(gamma%(A11*delta1));

double first = -0.50*dot(mean_piece, w%mean_piece) - 
               0.50*(1/sigma2_A)*(A11_trans*A11_trans);

/*Decision*/
double ratio = exp(first - second);   
double acc = 1;
if(ratio < R::runif(0, 1)){
  A11 = A11_old;
  acc = 0;
  }
acctot_A11_trans = acctot_A11_trans + 
                   acc;

return Rcpp::List::create(Rcpp::Named("A11") = A11,
                          Rcpp::Named("acctot_A11_trans") = acctot_A11_trans);

}



