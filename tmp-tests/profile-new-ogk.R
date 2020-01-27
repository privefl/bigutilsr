mat <- do.call("cbind", iris)
mat2 <- mat[rep(1:150, 5000), ]

Rcpp::sourceCpp("tmp-tests/ogk.cpp")

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

# profvis::profvis({test <- covRob_rcpp(mat2)})
