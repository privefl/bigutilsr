/******************************************************************************/
// Code originally written by Keurcien Luu
// See https://github.com/bcm-uga/pcadapt/blob/master/tmp-save/ogk.cpp
/******************************************************************************/

#include <Rcpp.h>
using namespace Rcpp;

inline double square(double x) { return x * x; }

/******************************************************************************/

// [[Rcpp::export]]
NumericVector scaleTau2_vector_rcpp(const NumericVector& x) {

  double c1 = 4.5;
  double c2 = 3.0;
  double c3 = 0.92471539217613152;
  double c4 = square(c2);

  int n = x.size();
  NumericVector x_dev(n);
  double medx = Rcpp::median(x);
  for (int i = 0; i < n; i++) {
    x_dev[i] = std::abs(x[i] - medx);
  }

  double sigma0 = Rcpp::median(x_dev);  // basically, MAD(x) (without constant)
  double sigma0_c1 = sigma0 * c1;
  double mu = 0;
  double w_tot = 0;
  for (int i = 0; i < n; i++) {
    double w_i = 1 - square(x_dev[i] / sigma0_c1);
    if (w_i > 0) {
      w_i *= w_i;
      mu += x[i] * w_i;
      w_tot += w_i;
    }
  }
  mu /= w_tot;

  double sum_rho = 0;
  for (int i = 0; i < n; i++) {
    double x_i_scaled = (x[i] - mu) / sigma0;
    sum_rho += std::min(square(x_i_scaled), c4);
  }
  sigma0 *= ::sqrt(sum_rho / (n * c3));

  return NumericVector::create(mu, sigma0);
}

/******************************************************************************/

// [[Rcpp::export]]
NumericVector dist_scaleTau2_matrix_rcpp(const NumericMatrix& Z) {

  int n = Z.nrow();
  int p = Z.ncol();
  NumericVector d(n);

  for (int j = 0; j < p; j++) {

    NumericVector col_j = Z.column(j);
    NumericVector mu_sigma0_j = scaleTau2_vector_rcpp(col_j);
    double mu_j     = mu_sigma0_j[0];
    double sigma0_j = mu_sigma0_j[1];

    for (int i = 0; i < n; i++) {
      d[i] += square((Z(i, j) - mu_j) / sigma0_j);
    }
  }

  return d;
}

/******************************************************************************/
