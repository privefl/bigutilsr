/******************************************************************************/
// Code originally written by Keurcien Luu
// See https://github.com/bcm-uga/pcadapt/blob/master/tmp-save/ogk.cpp
/******************************************************************************/

#include <Rcpp.h>
using namespace Rcpp;

inline double square(double x) { return x * x; }

/******************************************************************************/

// [[Rcpp::export]]
NumericVector& sum_in_temp(const NumericVector& x,
                           const NumericVector& y,
                           NumericVector& tmp_vec) {
  int n = x.size();
  for (int i = 0; i < n; i++) tmp_vec[i] = x[i] + y[i];

  return tmp_vec;
}

// [[Rcpp::export]]
NumericVector& sub_in_temp(const NumericVector& x,
                           const NumericVector& y,
                           NumericVector& tmp_vec) {
  int n = x.size();
  for (int i = 0; i < n; i++) tmp_vec[i] = x[i] - y[i];

  return tmp_vec;
}

/******************************************************************************/

double median_cpp(const NumericVector& x,
                  NumericVector& tmp_med) {

  std::copy(x.begin(), x.end(), tmp_med.begin());

  int n = tmp_med.size();
  bool is_even = (n % 2) == 0;
  int ind_mid = n / 2 - is_even;

  double * mid_ptr = tmp_med.begin() + ind_mid;
  std::nth_element(tmp_med.begin(), mid_ptr, tmp_med.end());
  double mid_val1 = tmp_med[ind_mid];
  if (!is_even) return mid_val1;

  std::nth_element(mid_ptr, mid_ptr + 1, tmp_med.end());
  double mid_val2 = mid_ptr[1];
  return (mid_val1 + mid_val2) / 2;
}

/******************************************************************************/

// [[Rcpp::export]]
NumericVector scaleTau2_vector_rcpp(const NumericVector& x,
                                    NumericVector& tmp_dev,
                                    NumericVector& tmp_med) {

  double c1 = 4.5;
  double c2 = 3.0;
  double c3 = 0.92471539217613152;
  double c4 = square(c2);

  int n = x.size();
  double medx = median_cpp(x, tmp_med);
  for (int i = 0; i < n; i++) {
    tmp_dev[i] = std::abs(x[i] - medx);
  }

  double sigma0 = median_cpp(tmp_dev, tmp_med);
  double sigma0_c1 = sigma0 * c1;
  double mu = 0;
  double w_tot = 0;
  for (int i = 0; i < n; i++) {
    double w_i = 1 - square(tmp_dev[i] / sigma0_c1);
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
  NumericVector d(n), tmp_dev(n), tmp_med(n);

  for (int j = 0; j < p; j++) {

    NumericVector col_j = Z.column(j);
    NumericVector mu_sigma0_j = scaleTau2_vector_rcpp(col_j, tmp_dev, tmp_med);
    double mu_j     = mu_sigma0_j[0];
    double sigma0_j = mu_sigma0_j[1];

    for (int i = 0; i < n; i++) {
      d[i] += square((Z(i, j) - mu_j) / sigma0_j);
    }
  }

  return d;
}

/******************************************************************************/
