// Code originally written by Keurcien Luu
// See https://github.com/bcm-uga/pcadapt/blob/master/tmp-save/ogk.cpp

#include <Rcpp.h>
using namespace Rcpp;

static const double c1 = 4.5;
static const double c2 = 3.0;
static const double c3 = 0.92471539217613152;

// [[Rcpp::export]]
NumericVector colMedian_rcpp(const NumericMatrix& x) {
  int p = x.ncol();
  NumericVector out(p);
  for (int j = 0; j < p; j++) {
    NumericVector tmp(x(_, j));
    out[j] = Rcpp::median(tmp);
  }
  return out;
}

List scaleTau2_matrix_rcpp(const NumericMatrix& x){
  int n = x.nrow();
  int p = x.ncol();
  NumericVector medx(p);
  NumericVector sigma0(p);
  medx = colMedian_rcpp(x);
  NumericVector xvec(n);
  NumericVector rho(n);
  NumericVector muvec(p);
  double mu, w_tot, sum_rho, tmp;
  for (int j = 0; j < p; j++) {
    mu = 0;
    w_tot = 0;
    sum_rho = 0;
    for (int i = 0; i < n; i++) {
      xvec[i] = std::abs(x(i, j) - medx[j]);
    }
    sigma0[j] = Rcpp::median(xvec);
    for (int i = 0; i < n; i++) {
      xvec[i] /= (sigma0[j] * c1);
      double w_i = 1 - (xvec[i] * xvec[i]);
      if (w_i > 0) {
        w_i *= w_i;
        mu += x(i, j) * w_i;
        w_tot += w_i;
      }
    }
    muvec[j] = mu / w_tot;
    for (int i = 0; i < n; i++) {
      tmp = ((x(i, j) - muvec[j]) / sigma0[j]) * ((x(i, j) - muvec[j]) / sigma0[j]);
      if (tmp > (c2 * c2)){
        rho[i] = c2 * c2;
      } else {
        rho[i] = tmp;
      }
      sum_rho += rho[i];
    }
    sigma0[j] *= ::sqrt(sum_rho / (n * c3));
  }
  return List::create(_["mu"] = muvec, _["sigma0"] = sigma0);
}

NumericVector scaleTau2_vector_rcpp(NumericVector &x) {
  int n = x.size();
  double medx;
  double sigma0;
  medx = median(x);
  NumericVector xvec(n);
  NumericVector rho(n);
  double mu = 0;
  double w_tot = 0;
  double sum_rho = 0;
  double tmp;
  for (int i = 0; i < n; i++) {
    xvec[i] = std::abs(x[i] - medx);
  }
  sigma0 = median(xvec);
  for (int i = 0; i < n; i++) {
    xvec[i] /= (sigma0 * c1);
    double w_i = 1 - (xvec[i] * xvec[i]);
    if (w_i > 0) {
      w_i *= w_i;
      mu += x[i] * w_i;
      w_tot += w_i;
    }
  }
  mu /= w_tot;
  for (int i = 0; i < n; i++) {
    tmp = ((x[i] - mu) / sigma0) * ((x[i] - mu) / sigma0);
    if (tmp > (c2 * c2)) {
      rho[i] = c2 * c2;
    } else {
      rho[i] = tmp;
    }
    sum_rho += rho[i];
  }
  sigma0 *= ::sqrt(sum_rho / (n * c3));
  return NumericVector::create(mu, sigma0);;
}

double covGK_rcpp(const NumericVector& x, const NumericVector& y){
  NumericVector sum_xy  = x + y;
  NumericVector diff_xy = x - y;
  NumericVector ms_sum  = scaleTau2_vector_rcpp(sum_xy);
  NumericVector ms_diff = scaleTau2_vector_rcpp(diff_xy);
  return(((ms_sum[1] * ms_sum[1]) - (ms_diff[1] * ms_diff[1])) / 4);
}

// [[Rcpp::export]]
List ogk_step_rcpp(const NumericMatrix& x) {

  int p = x.ncol();
  List ms = scaleTau2_matrix_rcpp(x);
  NumericVector sigma0 = ms[1];
  NumericMatrix U(p, p);
  for (int i = 0; i < p; i++) {
    for (int j = i; j < p; j++) {
      if (i == j) {
        U(i, j) = 1.0;
      } else {
        NumericVector tmp_i = x.column(i) / sigma0[i];
        NumericVector tmp_j = x.column(j) / sigma0[j];
        U(j, i) = U(i, j) = covGK_rcpp(tmp_i, tmp_j);
      }
    }
  }

  return List::create(_["U"] = U, _["sigma0"] = sigma0);
}

// [[Rcpp::export]]
List covRob_ogk_rcpp(const NumericMatrix& x, NumericMatrix& Z) {

  int n = x.nrow();
  int p = x.ncol();
  double df = (double) p;

  List musigma = scaleTau2_matrix_rcpp(Z);
  NumericVector mu = musigma[0];
  NumericVector sigma0 = musigma[1];
  NumericVector d(n);
  double medd;
  for (int j = 0; j < p; j++){
    for (int i = 0; i < n; i++){
      Z(i, j) -= mu[j];
      Z(i, j) /= sigma0[j];
      d[i] += Z(i, j) * Z(i, j);
    }
  }
  medd = Rcpp::median(d);
  double tmp = ::Rf_qchisq(0.5, df, 1, 0);
  double cdelta = medd / tmp;
  double beta = 0.9;
  double quantile = ::Rf_qchisq(beta, df, 1, 0);
  double qq = quantile * cdelta;
  double sum_wt = 0;

  NumericMatrix wcov(p, p);
  NumericVector wcenter(p);
  IntegerVector wt(n);

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

  for (int i = 0; i < p; i ++){
    for (int j = i; j < p; j++){
      for (int k = 0; k < n; k++){
        if (wt[k] == 1){
          wcov.at(i, j) += (x.at(k, i) - wcenter[i]) * (x.at(k, j) - wcenter[j]) / sum_wt;
        }
      }
      wcov.at(j, i) = wcov.at(i, j);
    }
  }

  return List::create(_["cov"] = wcov, _["center"] = wcenter);
}
