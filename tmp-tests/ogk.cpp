/******************************************************************************/
// Code originally written by Keurcien Luu
// See https://github.com/bcm-uga/pcadapt/blob/master/tmp-save/ogk.cpp
/******************************************************************************/

#include <Rcpp.h>
using namespace Rcpp;

static const double c1 = 4.5;
static const double c2 = 3.0;
static const double c3 = 0.92471539217613152;
static const double c4 = 9.0;

inline double square(double x) { return x * x; }

/******************************************************************************/

NumericVector scaleTau2_vector_rcpp(const NumericVector& x) {

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
List scaleTau2_matrix_rcpp(const NumericMatrix& x) {

  int p = x.ncol();
  NumericVector mu_vec(p), sigma0_vec(p);

  for (int j = 0; j < p; j++) {
    NumericVector col_j = x.column(j);
    NumericVector mu_sigma0_j = scaleTau2_vector_rcpp(col_j);
    mu_vec[j]     = mu_sigma0_j[0];
    sigma0_vec[j] = mu_sigma0_j[1];
  }

  return List::create(_["mu"] = mu_vec, _["sigma0"] = sigma0_vec);
}

/******************************************************************************/

// [[Rcpp::export]]
double covGK_rcpp(const NumericVector& x, const NumericVector& y){
  double sigma0_sum  = scaleTau2_vector_rcpp(x + y)[1];
  double sigma0_diff = scaleTau2_vector_rcpp(x - y)[1];
  return (square(sigma0_sum) - square(sigma0_diff)) / 4;
}

/******************************************************************************/

// [[Rcpp::export]]
double test_qchisq(double x, double df) {
  return ::Rf_qchisq(x, df, 1, 0);
}
