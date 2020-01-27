# mat <- bigutilsr::as_model_matrix(iris)
# mat <- bigutilsr::as_model_matrix(iris, intercept = TRUE)[, -1]
mat <- do.call("cbind", iris)
mat2 <- mat[rep(1:150, 500), ]

cov(mat2)

true <- robust::covRob(mat2, estim = "pairwiseGK")

if (!exists("tmp")) {
  download.file(
    "https://raw.githubusercontent.com/bcm-uga/pcadapt/master/tmp-save/ogk.cpp",
    tmp <- tempfile(fileext = ".cpp"))
}
Rcpp::sourceCpp(tmp)
test <- covRob_rcpp(mat2)
all.equal(true[c("cov", "center", "dist")], test[c("cov", "center", "dist")],
          check.attributes = FALSE)


Rcpp::sourceCpp("tmp-tests/ogk.cpp")
all.equal(test_qchisq(0.77, 2), qchisq(0.77, 2))
all.equal(test_qchisq(0.17, 8.77), qchisq(0.17, 8.77))

ogk_step_r <- function(X) {

  p <- ncol(X)
  U <- matrix(1, p, p)
  X.col.scaled <- list()
  sigma0 <- scaleTau2_matrix_rcpp(X)$sigma0

  for (j in seq_len(p)) {
    X.col.scaled[[j]] <- X[, j] / sigma0[j]
    for (k in seq_len(j - 1)) {
      U[j, k] <- U[k, j] <- covGK_rcpp(X.col.scaled[[j]], X.col.scaled[[k]]);
    }
  }

  eigvec <- eigen(U, symmetric = TRUE)$vectors
  X %*% sweep(eigvec, 1, sigma0, '/')
}

covRob_ogk_r <- function(X, beta = 0.9) {

  # First iteration
  V <- ogk_step_r(X)

  # Second iteration
  Z <- ogk_step_r(V)

  musigma0 <- scaleTau2_matrix_rcpp(Z);
  Z_scaled <- scale(Z, center = musigma0$mu, scale = musigma0$sigma0)
  d <- rowSums(Z_scaled ** 2)

  df <- ncol(X)
  cdelta <- median(d) / qchisq(0.5, df)
  quantile <- qchisq(beta, df)

  X.in <- X[d < (quantile * cdelta), ]
  wcenter <- colMeans(X.in)
  wcov <- crossprod(sweep(X.in, 2, wcenter, '-')) / nrow(X.in)

  list(cov = wcov, center = wcenter, dist = stats::mahalanobis(X, wcenter, wcov))
}
test2 <- covRob_ogk_r(mat2)
all.equal(true[c("cov", "center", "dist")], test2[c("cov", "center", "dist")],
          check.attributes = FALSE)

microbenchmark::microbenchmark(
  # ROB = robust::covRob(mat2, estim = "pairwiseGK"),
  KL = covRob_rcpp(mat2),
  FP = covRob_ogk_r(mat2),
  times = 15
)
# Unit: milliseconds
# expr      min       lq     mean   median       uq      max neval
#  ROB 676.8654 677.9064 683.2131 681.0307 684.9266 696.8045    10
#   KL 231.3390 233.4321 234.8431 234.5108 235.5153 240.4946    10
