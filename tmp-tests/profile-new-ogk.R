mat <- do.call("cbind", iris)
mat2 <- mat[rep(1:150, 5000), ]

Rcpp::sourceCpp("tmp-tests/ogk.cpp")

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
      sum <- col_j + col_k
      sum[1]
      sigma0_sum  <- get_sigma0(sum)
      diff <- col_j - col_k
      diff[1]
      sigma0_diff <- get_sigma0(diff)
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
profvis::profvis({test2 <- covRob_ogk_r(mat2)})

# profvis::profvis({test <- covRob_rcpp(mat2)})
