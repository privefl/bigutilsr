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
    eps <- sum((Z - Z_new) ** 2) / sum(sweep(Z_new, 2, colMeans(Z_new)) ** 2)
    # print(eps)
    if (eps < epsilon_min) break
    Z <- Z_new
  }
  proc
}

#### bed example ####
library(bigsnpr)
bedfile <- "~/Bureau/bigsnpr/tmp-data/1000G_phase3_common_norel.bed"
(X <- bed(bedfile))

N <- 1800; M <- ncol(X)
set.seed(1); ind <- sample(nrow(X), N)
K <- 20; K2 <- K; K3 <- K2
# K <- 14; K2 <- 14; K3 <- 20
obj.svd <- bed_autoSVD2(X, ind.row = ind, k = K3, ncores = nb_cores())
# saveRDS(obj.svd, "tmp-data/svd_1000G_K20.rds")
# obj.svd <- readRDS("tmp-data/svd_1000G_K20.rds")
U <- obj.svd$u
d <- obj.svd$d
V <- obj.svd$v

keep <- attr(obj.svd, "subset.col")
U2 <- big_apply(X, function(X, ind, ind.col, center, scale, V) {
  ind.row <- bigstatsr::rows_along(X)
  nb_nona <- bigsnpr:::bed_stats(X, ind.row, ind.col[ind])$nb_nona_row
  U <- apply(V[ind, , drop = FALSE], 2, function(v) {
    bigsnpr:::pMatVec4(X, ind.row, ind.col[ind], center[ind], scale[ind], v)
  })
  list(U, nb_nona)
}, ind = seq_along(keep), ncores = nb_cores(),
ind.col = keep,
center = obj.svd$center, scale = obj.svd$scale,
V = obj.svd$v)
U3 <- U2[[1]][[1]]
nb_nona <- U2[[1]][[2]]
for (u2 in U2[-1]) {
  U3 <- U3 + u2[[1]]
  nb_nona <- nb_nona + u2[[2]]
}
U3 <- sweep(U3, 1, length(keep) / nb_nona, '*')

plot(U3[, 5:6], col = (!rows_along(U3) %in% ind) + 1, pch = 20)

U1 <- U3[-ind, ]

PC.ref <- predict(obj.svd)[, 1:K]


U2 <- cbind(rbind(U, 0), 0)
U2[nrow(U2), ncol(U2)] <- 1
Q <- cbind(rbind(diag(d), 0), 0)
dim <- nrow(Q)

set.seed(2); ind2 <- sample(rows_along(X)[-ind], 200)
X2 <- bigsnpr:::read_bed_scaled(X, ind2, keep, center = obj.svd$center,
                                scale = obj.svd$scale)
proj <- matrix(0, length(ind2), K)
for (i in seq_along(ind2)) {
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

  proc <- procrustes_diffdim(PC.ref, PC.aug[-nrow(PC.aug), ])
  proj[i, ] <- apply_procrustes(PC.aug[nrow(PC.aug), , drop = FALSE], proc)[, 1:K]
}

plot(proj, U3[ind2, 1:K])

col <- 5:6
plot(PC.ref[, col])
points(U1[, col], col = "red", pch = 20)
points(proj[, col], col = "blue", pch = 20)
# points(proj2[, col], col = "pink", pch = 20)

k <- 12
plot(U3[ind2, k], proj[, k])
(mylm <- lm(proj[, k] ~ U3[ind2, k] + 0))
abline(mylm, col = "red")


coef2 <- sapply(1:K, function(j) {
  c(
    lm(proj[, j] ~ U3[ind2, j] + 0)$coefficients[[1]],
    #robust::lmRob(proj[, j] ~ U3[ind2, j] + 0)$coefficients[[1]],
    summary(lm(proj[, j] ~ U3[ind2, j] + 0))$r.squared,
    1 / lm(U3[ind2, j] ~ proj[, j] + 0)$coefficients[[1]],
    {
      pca.coef <- prcomp(cbind(proj[, j], U3[ind2, j]))$rot[, 1]
      pca.coef[1] / pca.coef[2]
    },
    #1 / robust::lmRob(U3[ind2, j] ~ proj[, j] + 0)$coefficients[[1]],
    NULL
  )
})
coef2

for (k in 1:10) {
  col <- 2 * k - 1:0 #10:11 + 7# c(8, 12)  #3:4 + 16
  plot(PC.ref[, col])
  points(U1[, col], col = "red", pch = 20)
  points(proj[, col], col = "blue", pch = 20)
  points(sweep(U1[, col], 2, coef2[1, col], '*'), col = "green", pch = 20)
}
col <- c(14, 16) #10:11 + 7# c(8, 12)  #3:4 + 16
plot(PC.ref[, col])
points(U1[, col], col = "red", pch = 20)
points(proj[, col], col = "blue", pch = 20)
points(sweep(U1[, col], 2, coef2[1, col], '*'), col = "green", pch = 20)

rbind(
  apply(PC.ref, 2, mad),
  apply(proj, 2, mad),
  apply(sweep(U1, 2, coef2[5, ], '*'), 2, mad),
  apply(sweep(U1, 2, coef2[4, ], '*'), 2, mad),
  apply(sweep(U1, 2, coef2[6, ], '*'), 2, mad)
)[, 1:K]
