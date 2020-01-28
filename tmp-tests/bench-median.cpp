#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double median_rcpp(const NumericVector& x) {
  return Rcpp::median(x);
}



/*** R
x <- rnorm(1e6)
microbenchmark::microbenchmark(
  median.default(x),
  median_rcpp(x)
)
x <- rnorm(1e8)
profvis::profvis({median.default(x)})
profvis::profvis({median_rcpp(x)})
*/
