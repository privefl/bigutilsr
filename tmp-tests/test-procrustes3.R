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
  if (ncol(X) == ncol(Y)) return(get_procrustes(Y, X))
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
K <- 10; K2 <- K; K3 <- 2 * K2; svd <- svd(X[ind, ], nu = K3, nv = K3)
U <- svd$u
d <- head(svd$d, K3)
V <- svd$v

U1 <- X[-ind, ] %*% V[, 1:K]
attr(bigutilsr::pca_adjust(U1, svd$d^2, M, N), "shrinkage")

PC.ref <- sweep(U, 2, d, '*')[, 1:K]


U2 <- cbind(rbind(U, 0), 0)
U2[nrow(U2), ncol(U2)] <- 1
Q <- cbind(rbind(diag(d), 0), 0)
dim <- nrow(Q)

X2 <- X[-ind, ]
proj <- U1
for (i in 1:nrow(proj)) {
  print(i)
  y <- X2[i, ]  # should be scaled
  L <- crossprod(V, y)
  H <- y - V %*% L
  H <- H / drop(sqrt(crossprod(H)))
  G <- crossprod(y, H)

  Q[, dim] <- c(L, G)
  R <- crossprod(Q)
  eig <- .Internal(La_rs(R, FALSE))  # eigen_real_symmetric_with_vectors
  PC.aug <- U2 %*% sweep(eig$vectors[, dim - 1:K2 + 1], 2,
                         sqrt(eig$values[dim - 1:K2 + 1]), '*')

  proc <- procrustes_diffdim(PC.ref, PC.aug[-301, ])
  proj[i, ] <- apply_procrustes(PC.aug[301, , drop = FALSE], proc)[, 1:K]
}

# proj2 <- t(apply(X[-ind, ], 1, function(x.new) {
#   X.aug <- rbind(x.new, X[ind, ])
#   svd.aug <- svd(X.aug, nu = K2, nv = 0)          # taking most of the time
#   PC.aug <- sweep(svd.aug$u, 2, svd.aug$d[1:K2], '*')
#   proc <- procrustes_diffdim(PC.ref, PC.aug[-1, ])
#   apply_procrustes(PC.aug[1, , drop = FALSE], proc)[, 1:K]
# }))

plot(proj, U1)
# plot(proj, proj2)

col <- 1:2
plot(PC.ref[, col])
points(U1[, col], col = "red", pch = 20)
points(proj[, col], col = "blue", pch = 20)
# points(proj2[, col], col = "pink", pch = 20)

# plot(U1[, 1], proj[, 1])
# (mylm <- lm(proj[, 1] ~ U1[, 1] + 0))
# abline(mylm, col = "red")
#
# plot(U1[, 2], proj[, 2])
# (mylm <- lm(proj[, 2] ~ U1[, 2] + 0))
# abline(mylm, col = "red")
#
# plot(U1[, 3], proj[, 3])
# (mylm <- lm(proj[, 3] ~ U1[, 3] + 0))
# abline(mylm, col = "red")
#
# plot(U1[, 4], proj[, 4])
# (mylm <- lm(proj[, 4] ~ U1[, 4] + 0))
# abline(mylm, col = "red")
#
# plot(U1[, 5], proj[, 5])
# (mylm <- lm(proj[, 5] ~ U1[, 5] + 0))
# abline(mylm, col = "red")

coef2 <- sapply(1:K, function(j) {
  c(
    lm(proj[, j] ~ U1[, j] + 0)$coefficients[[1]],
    robust::lmRob(proj[, j] ~ U1[, j] + 0)$coefficients[[1]],
    summary(lm(proj[, j] ~ U1[, j] + 0))$r.squared,
    # lm(proj2[, j] ~ U1[, j] + 0)$coefficients[[1]],
    1 / lm(U1[, j] ~ proj[, j] + 0)$coefficients[[1]],
    {
      pca.coef <- prcomp(cbind(proj[, j], U1[, j]))$rot[, 1]
      pca.coef[1] / pca.coef[2]
    },
    NULL
  )
})
coef2
# 1.1427442 1.3837123 1.8811833 2.2110356 2.2486710 2.3229277 2.3523555 2.443963
# hdpca::select.nspike(svd$d^2, M, N, n.spikes.max = 50)
(coef3 <- 1 / attr(bigutilsr::pca_adjust(U1, svd$d^2, M, N, n.spikes = 20), "shrinkage"))
# 1.146058 1.410460 2.153914 2.829381 2.983607 3.319876 3.550332

col <- 3:4 + 6
plot(PC.ref[, col])
points(U1[, col], col = "red", pch = 20)
# points(proj[, col], col = "blue", pch = 20)
points(sweep(U1[, col], 2, coef2[5, col], '*'), col = "green", pch = 20)
# points(sweep(U1[, col], 2, coef2[4, col], '*'), col = "pink", pch = 20)
points(sweep(U1[, col], 2, coef3[col], '*'), col = "pink", pch = 20)
points(proj[, col], col = "blue", pch = 20)

rbind(
  apply(PC.ref, 2, sd),
  apply(proj2, 2, sd),
  apply(proj, 2, sd),
  apply(sweep(U1, 2, coef2[5, ], '*'), 2, sd),
  apply(sweep(U1, 2, coef3, '*'), 2, sd)
)


# plot(proj[, col], sweep(U1[, col], 2, coef2[2, col], '*')[ind2, ]); abline(0, 1, col = "red")
