#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// cvxreg
List cvxreg(NumericVector x, NumericVector y);
RcppExport SEXP _Sshaped_cvxreg(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(cvxreg(x, y));
    return rcpp_result_gen;
END_RCPP
}
// sshapedreg
List sshapedreg(NumericVector x, NumericVector y);
RcppExport SEXP _Sshaped_sshapedreg(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(sshapedreg(x, y));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Sshaped_cvxreg", (DL_FUNC) &_Sshaped_cvxreg, 2},
    {"_Sshaped_sshapedreg", (DL_FUNC) &_Sshaped_sshapedreg, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_Sshaped(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
