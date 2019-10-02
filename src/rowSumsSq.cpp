#include <Rcpp.h>
using namespace Rcpp;

using std::size_t;

inline double square(double x) {
  return x * x;
}

// [[Rcpp::export]]
NumericVector rowSumsSq(const NumericMatrix& source) {

  size_t n = source.rows();
  size_t m = source.cols();
  size_t i, j;

  NumericVector res(n);

  for (j = 0; j < m; j++)
    for (i = 0; i < n; i++)
      res[i] += square(source(i, j));

  return res;
}
