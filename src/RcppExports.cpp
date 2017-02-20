// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// persistence_from_cover
Rcpp::NumericMatrix persistence_from_cover(Rcpp::S4 cover, int max_dim, Rcpp::String representation, Rcpp::String reduction);
RcppExport SEXP dcph_persistence_from_cover(SEXP coverSEXP, SEXP max_dimSEXP, SEXP representationSEXP, SEXP reductionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type cover(coverSEXP);
    Rcpp::traits::input_parameter< int >::type max_dim(max_dimSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type representation(representationSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type reduction(reductionSEXP);
    rcpp_result_gen = Rcpp::wrap(persistence_from_cover(cover, max_dim, representation, reduction));
    return rcpp_result_gen;
END_RCPP
}
