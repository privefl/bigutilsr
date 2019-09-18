X <- matrix(rnorm(12), 4)
Y <- 2 * X + rnorm(12, sd = 0.01) + 4

procrustes_centered <- function(Y.centered, X.centered,
                                X.norm = sum(X.centered ** 2)) {

  svd <- svd(crossprod(Y.centered, X.centered))
  rho <- sum(svd$d) / X.norm
  R <- tcrossprod(svd$v, svd$u)
  # tcrossprod(Y.centered, R) / rho
  rho * (X.centered %*% R)
}

procrustes <- function(Y, X) {

  X.cm <- colMeans(X)
  X.centered <- sweep(X, 2, X.cm)
  Y.centered <- sweep(Y, 2, colMeans(Y))
  # Y.new <- procrustes_centered(Y.centered, X.centered)
  # sweep(Y.new, 2, X.cm, '+')
  X.new <- procrustes_centered(Y.centered, X.centered)
  sweep(X.new, 2, colMeans(Y), '+')
}

# Y_new <- procrustes(Y, X); all.equal(Y_new, X)
X_new <- procrustes(Y, X); all.equal(X_new, Y)

procrustes_diffdim <- function(Y, X, n_iter_max = 10e3, epsilon_min = 1e-6) {

  Z <- matrix(0, nrow(Y), ncol(X) - ncol(Y))
  for (i in 1:n_iter_max) {
    W <- cbind(Y, Z)
    X_new <- procrustes(W, X)
    Z_new = X_new[, (ncol(Y) + 1):ncol(X), drop = FALSE]
    print(eps <- crossprod(Z - Z_new) / sum(sweep(Z_new, 2, colMeans(Z_new)) ** 2))
    if (eps < epsilon_min) break
    Z <- Z_new
  }
  X_new[, 1:ncol(Y)]
}
X2 <- cbind(X, rnorm(nrow(X)))
X_new <- procrustes_diffdim(Y, X2)
all.equal(X_new, Y)
