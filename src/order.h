#ifndef __ORDER_H_INCLUDED__   
#define __ORDER_H_INCLUDED__ 

#include <Rcpp.h>

Rcpp::IntegerVector order2(Rcpp::NumericVector x, Rcpp::NumericVector y);

Rcpp::IntegerVector order3(Rcpp::NumericVector x, Rcpp::NumericVector y, Rcpp::NumericVector z);

Rcpp::IntegerVector order4(Rcpp::NumericVector x1, Rcpp::NumericVector x2, Rcpp::NumericVector x3, Rcpp::NumericVector x4);

Rcpp::IntegerVector order5(Rcpp::NumericVector x1, Rcpp::NumericVector x2, Rcpp::NumericVector x3, Rcpp::NumericVector x4, Rcpp::NumericVector x5);

#endif  