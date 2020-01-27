mat <- do.call("cbind", iris)
mat2 <- mat[rep(1:150, 5000), ]

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
profvis::profvis({test2 <- covRob_ogk_r(mat2)})
