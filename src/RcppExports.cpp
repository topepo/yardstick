// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// yardstick_mcc_multiclass_impl
SEXP yardstick_mcc_multiclass_impl(SEXP C);
RcppExport SEXP _yardstick_yardstick_mcc_multiclass_impl(SEXP CSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type C(CSEXP);
    rcpp_result_gen = Rcpp::wrap(yardstick_mcc_multiclass_impl(C));
    return rcpp_result_gen;
END_RCPP
}
// yardstick_pr_curve_binary_impl
SEXP yardstick_pr_curve_binary_impl(SEXP truth, SEXP estimate, SEXP thresholds);
RcppExport SEXP _yardstick_yardstick_pr_curve_binary_impl(SEXP truthSEXP, SEXP estimateSEXP, SEXP thresholdsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type truth(truthSEXP);
    Rcpp::traits::input_parameter< SEXP >::type estimate(estimateSEXP);
    Rcpp::traits::input_parameter< SEXP >::type thresholds(thresholdsSEXP);
    rcpp_result_gen = Rcpp::wrap(yardstick_pr_curve_binary_impl(truth, estimate, thresholds));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_yardstick_yardstick_mcc_multiclass_impl", (DL_FUNC) &_yardstick_yardstick_mcc_multiclass_impl, 1},
    {"_yardstick_yardstick_pr_curve_binary_impl", (DL_FUNC) &_yardstick_yardstick_pr_curve_binary_impl, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_yardstick(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
