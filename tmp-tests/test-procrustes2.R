get_procrustes <- function(Y, X, X.norm = sum(X.centered ** 2)) {
  X_mean <- colMeans(X)
  X.centered <- sweep(X, 2, X_mean)
  svd <- svd(crossprod(Y, X.centered))
  rho <- sum(svd$d) / X.norm
  R <- tcrossprod(svd$v, svd$u)
  list(R = R, rho = rho, c = colMeans(Y) - crossprod(rho * X_mean, R))
}
apply_procrustes <- function(X, proc) {
  sweep(X %*% (proc$rho * proc$R), 2, proc$c, '+')
}
procrustes_diffdim <- function(Y, X, n_iter_max = 10e3, epsilon_min = 1e-6) {
  Z <- matrix(0, nrow(Y), ncol(X) - ncol(Y))
  for (i in 1:n_iter_max) {
    W <- cbind(Y, Z)
    proc <- get_procrustes(W, X)
    X_new <- apply_procrustes(X, proc)
    Z_new = X_new[, (ncol(Y) + 1):ncol(X), drop = FALSE]
    print(eps <- sum((Z - Z_new) ** 2) / sum(sweep(Z_new, 2, colMeans(Z_new)) ** 2))
    if (eps < epsilon_min) break
    Z <- Z_new
  }
  proc
}

#### First example ####
X <- readRDS(system.file("testdata", "three-pops.rds", package = "bigutilsr"))
pop <- rep(1:3, c(143, 167, 207))

N <- 300; M <- ncol(X)
ind <- sample(nrow(X), N)
svd <- svd(X[ind, ])
# plot(svd$d^2, log = "xy")
# hist(svd$d[svd$d < 80]^2, breaks = nclass.scottRob)

U1 <- X[-ind, ] %*% svd$v[, 1:4]
attr(bigutilsr::pca_adjust(U1, svd$d^2, M, N), "shrinkage")

PC.ref <- X[ind, ] %*% svd$v[, 1:4]
proj <- t(apply(X[-ind, ], 1, function(x.new) {
  X.aug <- rbind(x.new, X[ind, ])
  svd.aug <- svd(X.aug, nu = 10, nv = 0)
  PC.aug <- sweep(svd.aug$u, 2, svd.aug$d[1:10], '*')
  proc <- procrustes_diffdim(PC.ref, PC.aug[-1, ])
  apply_procrustes(PC.aug[1, , drop = FALSE], proc)
}))

rbind(
  apply(U1 / proj[, 1:4], 2, median),
  attr(bigutilsr::pca_adjust(U1, svd$d^2, M, N, n.spikes = 4), "shrinkage")
)

col <- 3:4
plot(PC.ref[, col])
points(U1[, col], col = "red", pch = 20)
points(proj[, col], col = "blue", pch = 20)
