#include "RcppArmadillo.h"
#include "CWVS.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double neg_two_loglike_update(arma::vec y,
                              arma::mat x,
                              arma::mat z, 
                              arma::vec off_set,
                              int likelihood_indicator,
                              int r,
                              double sigma2_epsilon,
                              arma::vec beta,
                              arma::vec gamma,
                              double A11,
                              arma::vec delta1){

int n = y.size();
arma::vec dens(n); dens.fill(0.00);

arma::vec mu = off_set +
               x*beta + 
               z*(gamma%(A11*delta1));

if(likelihood_indicator == 0){
  
  arma::vec probs = exp(mu)/(1.00 + exp(mu));

  for(int j = 0; j < n; ++j){
     dens(j) = R::dbinom(y(j),
                         1,
                         probs(j),
                         TRUE);
     }

  }

if(likelihood_indicator == 1){
  for(int j = 0; j < n; ++j){
     dens(j) = R::dnorm(y(j),
                        mu(j),
                        sqrt(sigma2_epsilon),
                        TRUE);
     }
  }

if(likelihood_indicator == 2){
  
  arma::vec probs = exp(mu)/(1.00 + exp(mu));
  for(int j = 0; j < n; ++j){
     dens(j) = R::dnbinom(y(j), 
                          r, 
                          (1.00 - probs(j)),        
                          TRUE);
     }
  
  }

double neg_two_loglike = -2.00*sum(dens);

return neg_two_loglike;

}

























































