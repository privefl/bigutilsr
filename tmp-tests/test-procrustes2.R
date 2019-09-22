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
X <- scale(X)
pop <- rep(1:3, c(143, 167, 207))

N <- 300; M <- ncol(X)
ind <- sample(nrow(X), N)
svd <- svd(X[ind, ])
# plot(svd$d^2, log = "xy")
# hist(svd$d[svd$d < 80]^2, breaks = nclass.scottRob)

U1 <- X[-ind, ] %*% svd$v[, 1:10]
attr(bigutilsr::pca_adjust(U1, svd$d^2, M, N), "shrinkage")

PC.ref <- X[ind, ] %*% svd$v[, 1:5]

ind2 <- sample(217, 50)
proj <- t(apply(X[-ind, ][ind2, ], 1, function(x.new) {
  X.aug <- rbind(x.new, X[ind, ])
  svd.aug <- svd(X.aug, nu = 10, nv = 0)          # taking most of the time
  PC.aug <- sweep(svd.aug$u, 2, svd.aug$d[1:10], '*')
  proc <- procrustes_diffdim(PC.ref, PC.aug[-1, ])
  apply_procrustes(PC.aug[1, , drop = FALSE], proc)
}))


(coef <- apply(U1[ind2, ] / proj, 2, median))
(attr(bigutilsr::pca_adjust(U1, svd$d^2, M, N, n.spikes = 5), "shrinkage"))

col <- 4:5
plot(PC.ref[, col])
points(U1[, col], col = "red", pch = 20)
points(proj[, col], col = "blue", pch = 20)
points(sweep(U1[, col], 2, coef[col], '/'), col = "green", pch = 20)


plot(U1[ind2, 1], proj[, 1])
(mylm <- lm(proj[, 1] ~ U1[ind2, 1] + 0))
abline(mylm, col = "red")

plot(U1[ind2, 2], proj[, 2])
(mylm <- lm(proj[, 2] ~ U1[ind2, 2] + 0))
abline(mylm, col = "red")

plot(U1[ind2, 3], proj[, 3])
(mylm <- lm(proj[, 3] ~ U1[ind2, 3] + 0))
abline(mylm, col = "red")

plot(U1[ind2, 4], proj[, 4])
(mylm <- lm(proj[, 4] ~ U1[ind2, 4] + 0))
abline(mylm, col = "red")

plot(U1[ind2, 5], proj[, 5])
(mylm <- lm(proj[, 5] ~ U1[ind2, 5] + 0))
abline(mylm, col = "red")

coef2 <- sapply(1:10, function(j) {
  c(
    lm(proj[, j] ~ U1[ind2, j] + 0)$coefficients[[1]],
    robust::lmRob(proj[, j] ~ U1[ind2, j] + 0)$coefficients[[1]],
    summary(lm(proj[, j] ~ U1[ind2, j] + 0))$r.squared
  )
})
coef2

col <- 1:2
plot(PC.ref[, col])
points(U1[, col], col = "red", pch = 20)
# points(proj[, col], col = "blue", pch = 20)
points(sweep(U1[, col], 2, coef2[2, col], '*'), col = "green", pch = 20)

# plot(proj[, col], sweep(U1[, col], 2, coef2[2, col], '*')[ind2, ]); abline(0, 1, col = "red")
