#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List find_nn(const NumericMatrix& X, const IntegerVector& ind) {

  int N = X.nrow();
  int M = X.ncol();
  int K = ind.size();
  IntegerVector I(K);
  NumericVector D(K);

  for (int k = 0; k < K; k++) {
    int i_k = ind[k] - 1;
    int i_min = -1;
    double d_min = R_PosInf;
    for (int i = 0; i < N; i++) {
      if (i == i_k) continue;
      double d_max = 0;
      for (int j = 0; j < M; j++) {
        double d = std::abs(X(i, j) - X(i_k, j));
        if (d > d_max) {
          d_max = d;
          if (d_max >= d_min) break;
        }
      }
      if (d_max < d_min) {
        i_min = i;
        d_min = d_max;
      }
    }
    I[k] = i_min + 1;
    D[k] = d_min;
  }

  return List::create(_["id"] = I, _["dist"] = D);
}


/*** R
X <- matrix(rnorm(2000), ncol = 5)
find_nn(X, 3:4)
nabor::knn(X, X[3:4, ], k = 2)
microbenchmark::microbenchmark(
  find_nn(X, 3:4),
  nabor::knn(X, X[3:4, ], k = 2)
)
*/
