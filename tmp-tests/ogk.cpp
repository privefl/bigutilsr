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

double covGK_rcpp(const NumericVector& x, const NumericVector& y){
  double sigma0_sum  = scaleTau2_vector_rcpp(x + y)[1];
  double sigma0_diff = scaleTau2_vector_rcpp(x - y)[1];
  return (square(sigma0_sum) - square(sigma0_diff)) / 4;
}

/******************************************************************************/

// [[Rcpp::export]]
NumericMatrix ogk_step_rcpp(const NumericMatrix& x_scaled) {

  int p = x_scaled.ncol();
  NumericMatrix U(p, p);

  for (int i = 0; i < p; i++) {
    NumericVector col_i = x_scaled.column(i);
    for (int j = 0; j < i; j++) {
      NumericVector col_j = x_scaled.column(j);
      U(j, i) = U(i, j) = covGK_rcpp(col_i, col_j);
    }
    U(i, i) = 1.0;
  }

  return U;
}

/******************************************************************************/

// [[Rcpp::export]]
List covRob_ogk_rcpp(const NumericMatrix& x, NumericMatrix& Z) {

  int n = x.nrow();
  int p = x.ncol();
  double df = (double) p;

  List musigma = scaleTau2_matrix_rcpp(Z);
  NumericVector mu     = musigma[0];
  NumericVector sigma0 = musigma[1];
  NumericVector d(n);

  for (int j = 0; j < p; j++) {
    for (int i = 0; i < n; i++) {
      Z(i, j) -= mu[j];
      Z(i, j) /= sigma0[j];
      d[i] += square(Z(i, j));
    }
  }

  double cdelta = Rcpp::median(d) / ::Rf_qchisq(0.5, df, 1, 0);
  double beta = 0.9;
  double quantile = ::Rf_qchisq(beta, df, 1, 0);
  double qq = quantile * cdelta;

  NumericVector wcenter(p);
  IntegerVector wt(n);
  double sum_wt = 0;

  for (int i = 0; i < n; i++) {
    if (d[i] < qq){
      wt[i] = 1;
      sum_wt += 1;
      for (int j = 0; j < p; j++) {
        wcenter[j] += x.at(i, j);
      }
    }
  }
  for (int j = 0; j < p; j++) {
    wcenter[j] /= sum_wt;
  }

  NumericMatrix wcov(p, p);

  for (int i = 0; i < p; i ++){
    for (int j = i; j < p; j++){
      for (int k = 0; k < n; k++){
        if (wt[k] == 1){
          wcov.at(i, j) += (x.at(k, i) - wcenter[i]) * (x.at(k, j) - wcenter[j]);
        }
      }
      wcov.at(i, j) /= sum_wt;
      wcov.at(j, i) = wcov.at(i, j);
    }
  }

  return List::create(_["cov"] = wcov, _["center"] = wcenter);
}

/******************************************************************************/
