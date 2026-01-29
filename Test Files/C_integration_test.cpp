#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
int addC(int x, int y) {
  int sum = x + y;
  return sum;
}