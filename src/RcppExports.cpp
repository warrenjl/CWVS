// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// A11_update
Rcpp::List A11_update(double A11_old, arma::mat x, arma::mat z, arma::vec w, arma::vec gamma_l, arma::vec beta, arma::vec gamma, arma::vec delta1, double sigma2_A, double metrop_var_A11_trans, int acctot_A11_trans);
RcppExport SEXP _CWVS_A11_update(SEXP A11_oldSEXP, SEXP xSEXP, SEXP zSEXP, SEXP wSEXP, SEXP gamma_lSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP delta1SEXP, SEXP sigma2_ASEXP, SEXP metrop_var_A11_transSEXP, SEXP acctot_A11_transSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type A11_old(A11_oldSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma_l(gamma_lSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type delta1(delta1SEXP);
    Rcpp::traits::input_parameter< double >::type sigma2_A(sigma2_ASEXP);
    Rcpp::traits::input_parameter< double >::type metrop_var_A11_trans(metrop_var_A11_transSEXP);
    Rcpp::traits::input_parameter< int >::type acctot_A11_trans(acctot_A11_transSEXP);
    rcpp_result_gen = Rcpp::wrap(A11_update(A11_old, x, z, w, gamma_l, beta, gamma, delta1, sigma2_A, metrop_var_A11_trans, acctot_A11_trans));
    return rcpp_result_gen;
END_RCPP
}
// A21_update
double A21_update(arma::vec gamma_star, arma::vec delta1, double A22, arma::vec delta2, double sigma2_A);
RcppExport SEXP _CWVS_A21_update(SEXP gamma_starSEXP, SEXP delta1SEXP, SEXP A22SEXP, SEXP delta2SEXP, SEXP sigma2_ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type gamma_star(gamma_starSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type delta1(delta1SEXP);
    Rcpp::traits::input_parameter< double >::type A22(A22SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type delta2(delta2SEXP);
    Rcpp::traits::input_parameter< double >::type sigma2_A(sigma2_ASEXP);
    rcpp_result_gen = Rcpp::wrap(A21_update(gamma_star, delta1, A22, delta2, sigma2_A));
    return rcpp_result_gen;
END_RCPP
}
// A22_update
Rcpp::List A22_update(double A22_old, arma::vec gamma_star, arma::vec delta1, double A21_old, arma::vec delta2, double sigma2_A, double metrop_var_A22_trans, int acctot_A22_trans);
RcppExport SEXP _CWVS_A22_update(SEXP A22_oldSEXP, SEXP gamma_starSEXP, SEXP delta1SEXP, SEXP A21_oldSEXP, SEXP delta2SEXP, SEXP sigma2_ASEXP, SEXP metrop_var_A22_transSEXP, SEXP acctot_A22_transSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type A22_old(A22_oldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma_star(gamma_starSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type delta1(delta1SEXP);
    Rcpp::traits::input_parameter< double >::type A21_old(A21_oldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type delta2(delta2SEXP);
    Rcpp::traits::input_parameter< double >::type sigma2_A(sigma2_ASEXP);
    Rcpp::traits::input_parameter< double >::type metrop_var_A22_trans(metrop_var_A22_transSEXP);
    Rcpp::traits::input_parameter< int >::type acctot_A22_trans(acctot_A22_transSEXP);
    rcpp_result_gen = Rcpp::wrap(A22_update(A22_old, gamma_star, delta1, A21_old, delta2, sigma2_A, metrop_var_A22_trans, acctot_A22_trans));
    return rcpp_result_gen;
END_RCPP
}
// CWVS
Rcpp::List CWVS(int mcmc_samples, arma::vec y, arma::mat x, arma::mat z, double metrop_var_phi1_trans, double metrop_var_phi2_trans, double metrop_var_A11_trans, double metrop_var_A22_trans, Rcpp::Nullable<double> sigma2_beta_prior, Rcpp::Nullable<double> alpha_phi1_prior, Rcpp::Nullable<double> beta_phi1_prior, Rcpp::Nullable<double> alpha_phi2_prior, Rcpp::Nullable<double> beta_phi2_prior, Rcpp::Nullable<double> sigma2_A_prior, Rcpp::Nullable<Rcpp::NumericVector> beta_init, Rcpp::Nullable<Rcpp::NumericVector> gamma_init, Rcpp::Nullable<Rcpp::NumericVector> delta1_init, Rcpp::Nullable<Rcpp::NumericVector> delta2_init, Rcpp::Nullable<double> phi1_init, Rcpp::Nullable<double> phi2_init, Rcpp::Nullable<double> A11_init, Rcpp::Nullable<double> A22_init, Rcpp::Nullable<double> A21_init);
RcppExport SEXP _CWVS_CWVS(SEXP mcmc_samplesSEXP, SEXP ySEXP, SEXP xSEXP, SEXP zSEXP, SEXP metrop_var_phi1_transSEXP, SEXP metrop_var_phi2_transSEXP, SEXP metrop_var_A11_transSEXP, SEXP metrop_var_A22_transSEXP, SEXP sigma2_beta_priorSEXP, SEXP alpha_phi1_priorSEXP, SEXP beta_phi1_priorSEXP, SEXP alpha_phi2_priorSEXP, SEXP beta_phi2_priorSEXP, SEXP sigma2_A_priorSEXP, SEXP beta_initSEXP, SEXP gamma_initSEXP, SEXP delta1_initSEXP, SEXP delta2_initSEXP, SEXP phi1_initSEXP, SEXP phi2_initSEXP, SEXP A11_initSEXP, SEXP A22_initSEXP, SEXP A21_initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type mcmc_samples(mcmc_samplesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type metrop_var_phi1_trans(metrop_var_phi1_transSEXP);
    Rcpp::traits::input_parameter< double >::type metrop_var_phi2_trans(metrop_var_phi2_transSEXP);
    Rcpp::traits::input_parameter< double >::type metrop_var_A11_trans(metrop_var_A11_transSEXP);
    Rcpp::traits::input_parameter< double >::type metrop_var_A22_trans(metrop_var_A22_transSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type sigma2_beta_prior(sigma2_beta_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type alpha_phi1_prior(alpha_phi1_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type beta_phi1_prior(beta_phi1_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type alpha_phi2_prior(alpha_phi2_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type beta_phi2_prior(beta_phi2_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type sigma2_A_prior(sigma2_A_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type beta_init(beta_initSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type gamma_init(gamma_initSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type delta1_init(delta1_initSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type delta2_init(delta2_initSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type phi1_init(phi1_initSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type phi2_init(phi2_initSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type A11_init(A11_initSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type A22_init(A22_initSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type A21_init(A21_initSEXP);
    rcpp_result_gen = Rcpp::wrap(CWVS(mcmc_samples, y, x, z, metrop_var_phi1_trans, metrop_var_phi2_trans, metrop_var_A11_trans, metrop_var_A22_trans, sigma2_beta_prior, alpha_phi1_prior, beta_phi1_prior, alpha_phi2_prior, beta_phi2_prior, sigma2_A_prior, beta_init, gamma_init, delta1_init, delta2_init, phi1_init, phi2_init, A11_init, A22_init, A21_init));
    return rcpp_result_gen;
END_RCPP
}
// beta_update
arma::vec beta_update(arma::mat x, arma::mat z, double sigma2_beta, arma::vec w, arma::vec gamma_l, arma::vec gamma_old, double A11_old, arma::vec delta1_old);
RcppExport SEXP _CWVS_beta_update(SEXP xSEXP, SEXP zSEXP, SEXP sigma2_betaSEXP, SEXP wSEXP, SEXP gamma_lSEXP, SEXP gamma_oldSEXP, SEXP A11_oldSEXP, SEXP delta1_oldSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2_beta(sigma2_betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma_l(gamma_lSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma_old(gamma_oldSEXP);
    Rcpp::traits::input_parameter< double >::type A11_old(A11_oldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type delta1_old(delta1_oldSEXP);
    rcpp_result_gen = Rcpp::wrap(beta_update(x, z, sigma2_beta, w, gamma_l, gamma_old, A11_old, delta1_old));
    return rcpp_result_gen;
END_RCPP
}
// delta1_update
arma::vec delta1_update(arma::mat x, arma::mat z, arma::vec w, arma::vec gamma_l, arma::vec beta, arma::vec gamma, arma::vec gamma_star, double A11_old, double A21_old, double A22_old, arma::vec delta2_old, arma::mat corr_inv1);
RcppExport SEXP _CWVS_delta1_update(SEXP xSEXP, SEXP zSEXP, SEXP wSEXP, SEXP gamma_lSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP gamma_starSEXP, SEXP A11_oldSEXP, SEXP A21_oldSEXP, SEXP A22_oldSEXP, SEXP delta2_oldSEXP, SEXP corr_inv1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma_l(gamma_lSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma_star(gamma_starSEXP);
    Rcpp::traits::input_parameter< double >::type A11_old(A11_oldSEXP);
    Rcpp::traits::input_parameter< double >::type A21_old(A21_oldSEXP);
    Rcpp::traits::input_parameter< double >::type A22_old(A22_oldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type delta2_old(delta2_oldSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type corr_inv1(corr_inv1SEXP);
    rcpp_result_gen = Rcpp::wrap(delta1_update(x, z, w, gamma_l, beta, gamma, gamma_star, A11_old, A21_old, A22_old, delta2_old, corr_inv1));
    return rcpp_result_gen;
END_RCPP
}
// delta2_update
arma::vec delta2_update(arma::vec gamma_star, arma::vec delta1, double A21_old, double A22_old, arma::mat corr_inv2);
RcppExport SEXP _CWVS_delta2_update(SEXP gamma_starSEXP, SEXP delta1SEXP, SEXP A21_oldSEXP, SEXP A22_oldSEXP, SEXP corr_inv2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type gamma_star(gamma_starSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type delta1(delta1SEXP);
    Rcpp::traits::input_parameter< double >::type A21_old(A21_oldSEXP);
    Rcpp::traits::input_parameter< double >::type A22_old(A22_oldSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type corr_inv2(corr_inv2SEXP);
    rcpp_result_gen = Rcpp::wrap(delta2_update(gamma_star, delta1, A21_old, A22_old, corr_inv2));
    return rcpp_result_gen;
END_RCPP
}
// gamma_star_update
arma::vec gamma_star_update(arma::vec gamma, arma::vec delta1_old, double A21_old, double A22_old, arma::vec delta2_old);
RcppExport SEXP _CWVS_gamma_star_update(SEXP gammaSEXP, SEXP delta1_oldSEXP, SEXP A21_oldSEXP, SEXP A22_oldSEXP, SEXP delta2_oldSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type delta1_old(delta1_oldSEXP);
    Rcpp::traits::input_parameter< double >::type A21_old(A21_oldSEXP);
    Rcpp::traits::input_parameter< double >::type A22_old(A22_oldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type delta2_old(delta2_oldSEXP);
    rcpp_result_gen = Rcpp::wrap(gamma_star_update(gamma, delta1_old, A21_old, A22_old, delta2_old));
    return rcpp_result_gen;
END_RCPP
}
// gamma_update
arma::vec gamma_update(arma::mat x, arma::mat z, arma::vec w, arma::vec gamma_l, arma::vec beta, arma::vec gamma_old, double A11_old, arma::vec delta1_old, double A21_old, double A22_old, arma::vec delta2_old);
RcppExport SEXP _CWVS_gamma_update(SEXP xSEXP, SEXP zSEXP, SEXP wSEXP, SEXP gamma_lSEXP, SEXP betaSEXP, SEXP gamma_oldSEXP, SEXP A11_oldSEXP, SEXP delta1_oldSEXP, SEXP A21_oldSEXP, SEXP A22_oldSEXP, SEXP delta2_oldSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma_l(gamma_lSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma_old(gamma_oldSEXP);
    Rcpp::traits::input_parameter< double >::type A11_old(A11_oldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type delta1_old(delta1_oldSEXP);
    Rcpp::traits::input_parameter< double >::type A21_old(A21_oldSEXP);
    Rcpp::traits::input_parameter< double >::type A22_old(A22_oldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type delta2_old(delta2_oldSEXP);
    rcpp_result_gen = Rcpp::wrap(gamma_update(x, z, w, gamma_l, beta, gamma_old, A11_old, delta1_old, A21_old, A22_old, delta2_old));
    return rcpp_result_gen;
END_RCPP
}
// neg_two_loglike_update
double neg_two_loglike_update(arma::vec y, arma::mat x, arma::mat z, arma::vec beta, arma::vec gamma, double A11, arma::vec delta1);
RcppExport SEXP _CWVS_neg_two_loglike_update(SEXP ySEXP, SEXP xSEXP, SEXP zSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP A11SEXP, SEXP delta1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type A11(A11SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type delta1(delta1SEXP);
    rcpp_result_gen = Rcpp::wrap(neg_two_loglike_update(y, x, z, beta, gamma, A11, delta1));
    return rcpp_result_gen;
END_RCPP
}
// phi_update
Rcpp::List phi_update(double phi_old, arma::vec delta, Rcpp::List temporal_corr_info, double alpha_phi, double beta_phi, double metrop_var_phi_trans, int acctot_phi_trans);
RcppExport SEXP _CWVS_phi_update(SEXP phi_oldSEXP, SEXP deltaSEXP, SEXP temporal_corr_infoSEXP, SEXP alpha_phiSEXP, SEXP beta_phiSEXP, SEXP metrop_var_phi_transSEXP, SEXP acctot_phi_transSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type phi_old(phi_oldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type temporal_corr_info(temporal_corr_infoSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_phi(alpha_phiSEXP);
    Rcpp::traits::input_parameter< double >::type beta_phi(beta_phiSEXP);
    Rcpp::traits::input_parameter< double >::type metrop_var_phi_trans(metrop_var_phi_transSEXP);
    Rcpp::traits::input_parameter< int >::type acctot_phi_trans(acctot_phi_transSEXP);
    rcpp_result_gen = Rcpp::wrap(phi_update(phi_old, delta, temporal_corr_info, alpha_phi, beta_phi, metrop_var_phi_trans, acctot_phi_trans));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_pgdraw
arma::vec rcpp_pgdraw(double b, arma::vec c);
RcppExport SEXP _CWVS_rcpp_pgdraw(SEXP bSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_pgdraw(b, c));
    return rcpp_result_gen;
END_RCPP
}
// temporal_corr_fun
Rcpp::List temporal_corr_fun(int p_z, double phi);
RcppExport SEXP _CWVS_temporal_corr_fun(SEXP p_zSEXP, SEXP phiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type p_z(p_zSEXP);
    Rcpp::traits::input_parameter< double >::type phi(phiSEXP);
    rcpp_result_gen = Rcpp::wrap(temporal_corr_fun(p_z, phi));
    return rcpp_result_gen;
END_RCPP
}
// w_update
Rcpp::List w_update(arma::vec y, arma::mat x, arma::mat z, arma::vec beta_old, arma::vec gamma_old, double A11_old, arma::vec delta1_old);
RcppExport SEXP _CWVS_w_update(SEXP ySEXP, SEXP xSEXP, SEXP zSEXP, SEXP beta_oldSEXP, SEXP gamma_oldSEXP, SEXP A11_oldSEXP, SEXP delta1_oldSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta_old(beta_oldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma_old(gamma_oldSEXP);
    Rcpp::traits::input_parameter< double >::type A11_old(A11_oldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type delta1_old(delta1_oldSEXP);
    rcpp_result_gen = Rcpp::wrap(w_update(y, x, z, beta_old, gamma_old, A11_old, delta1_old));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CWVS_A11_update", (DL_FUNC) &_CWVS_A11_update, 11},
    {"_CWVS_A21_update", (DL_FUNC) &_CWVS_A21_update, 5},
    {"_CWVS_A22_update", (DL_FUNC) &_CWVS_A22_update, 8},
    {"_CWVS_CWVS", (DL_FUNC) &_CWVS_CWVS, 23},
    {"_CWVS_beta_update", (DL_FUNC) &_CWVS_beta_update, 8},
    {"_CWVS_delta1_update", (DL_FUNC) &_CWVS_delta1_update, 12},
    {"_CWVS_delta2_update", (DL_FUNC) &_CWVS_delta2_update, 5},
    {"_CWVS_gamma_star_update", (DL_FUNC) &_CWVS_gamma_star_update, 5},
    {"_CWVS_gamma_update", (DL_FUNC) &_CWVS_gamma_update, 11},
    {"_CWVS_neg_two_loglike_update", (DL_FUNC) &_CWVS_neg_two_loglike_update, 7},
    {"_CWVS_phi_update", (DL_FUNC) &_CWVS_phi_update, 7},
    {"_CWVS_rcpp_pgdraw", (DL_FUNC) &_CWVS_rcpp_pgdraw, 2},
    {"_CWVS_temporal_corr_fun", (DL_FUNC) &_CWVS_temporal_corr_fun, 2},
    {"_CWVS_w_update", (DL_FUNC) &_CWVS_w_update, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_CWVS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
