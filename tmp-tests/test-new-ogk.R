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
  tmp <- ogk_step_rcpp(X)
  eigvec <- eigen(tmp$U, symmetric = TRUE)$vectors
  X %*% sweep(eigvec, 1, tmp$sigma0, '/')
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
  LUU = covRob_rcpp(mat2),
  LUU2 = covRob_ogk_r(mat2),
  times = 10
)
# Unit: milliseconds
# expr      min       lq     mean   median       uq      max neval
#  ROB 676.8654 677.9064 683.2131 681.0307 684.9266 696.8045    10
#  LUU 231.3390 233.4321 234.8431 234.5108 235.5153 240.4946    10

# microbenchmark::microbenchmark(
#   MED = apply(mat2, 2, median),
#   MED2 = colMedian_rcpp(mat2),
#   MED3 = colMedian_rcpp2(mat2),
#   times = 20
# )
# Unit: milliseconds
# expr    min      lq      mean  median      uq     max neval
#  MED 8.9783 9.23055 11.982365 9.56475 16.2671 18.3489    20
# MED2 4.2670 4.39190  4.844825 4.46860  4.6000 11.7844    20
# MED3 4.2349 4.31610  5.104240 4.43000  4.6205 10.8087    20
