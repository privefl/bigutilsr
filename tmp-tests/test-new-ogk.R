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

get_sigma0 <- function(x) {
  scaleTau2_vector_rcpp(x)[2]
}

ogk_step_r <- function(X) {

  p <- ncol(X)
  U <- matrix(1, p, p)
  X.col.scaled <- vector("list", p)

  for (j in seq_len(p)) {
    X_j <- X[, j]
    sigma0_j <- get_sigma0(X_j)
    col_j <- X.col.scaled[[j]] <- X_j / sigma0_j
    for (k in seq_len(j - 1)) {
      col_k <- X.col.scaled[[k]]
      sigma0_sum  <- get_sigma0(col_j + col_k)
      sigma0_diff <- get_sigma0(col_j - col_k)
      U[j, k] <- U[k, j] <- (sigma0_sum ** 2 - sigma0_diff ** 2) / 4
    }
  }

  do.call("cbind", X.col.scaled) %*% eigen(U, symmetric = TRUE)$vectors
}

covRob_ogk_r <- function(X, beta = 0.9) {

  # First iteration
  V <- ogk_step_r(X)

  # Second iteration
  Z <- ogk_step_r(V)

  d <- dist_scaleTau2_matrix_rcpp(Z)

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
# expr       min        lq      mean    median        uq       max neval
#  ROB 1308.2220 1458.1762 1604.0048 1548.4730 1731.3024 2001.8148    15
#   KL  435.0238  500.1638  608.0526  528.7331  654.5196 1052.8057    15
#   FP  268.5768  288.3206  358.0729  343.0261  369.8566  607.4272    15
