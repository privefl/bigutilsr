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

ogk_step_r <- function(X) {
  sigma0 <- scaleTau2_matrix_rcpp(X)$sigma0
  U <- ogk_step_rcpp(sweep(X, 2, sigma0, '/'))
  eigvec <- eigen(U, symmetric = TRUE)$vectors
  X %*% sweep(eigvec, 1, sigma0, '/')
}

covRob_ogk_r <- function(X) {

  # First iteration
  V <- ogk_step_r(X)

  # Second iteration
  Z <- ogk_step_r(V)
  res <- covRob_ogk_rcpp(X, Z)

  # Distance computation
  res$dist <- stats::mahalanobis(X, res$center, res$cov)

  res
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
