// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// addSimplices
IntegerVector addSimplices(IntegerVector x, IntegerVector y);
RcppExport SEXP hclustPH_addSimplices(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type y(ySEXP);
    __result = Rcpp::wrap(addSimplices(x, y));
    return __result;
END_RCPP
}
