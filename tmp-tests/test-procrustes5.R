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

#### bed example ####
library(bigsnpr)
bedfile <- download_1000G("tmp-data")
(X <- bed(bedfile))

N <- 1800; M <- ncol(X)
set.seed(1); ind <- sample(nrow(X), N)
K <- 30
obj.svd <- bed_autoSVD2(X, ind.row = ind, k = K, ncores = nb_cores())
U <- obj.svd$u
d <- obj.svd$d
V <- obj.svd$v

keep <- attr(obj.svd, "subset.col")
{
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
}

plot(U3[, 5:6], col = (!rows_along(U3) %in% ind) + 1, pch = 20)

U1 <- U3[-ind, ]

PC.ref <- predict(obj.svd)


U4 <- cbind(rbind(U, 0), 0)
U4[nrow(U4), ncol(U4)] <- 1
Q <- cbind(rbind(diag(d), 0), 0)
dim <- nrow(Q)

set.seed(2); ind2 <- rows_along(X)[-ind] #sample(rows_along(X)[-ind], 200)
X2 <- bigsnpr:::read_bed_scaled(X, ind2, keep,
                                center = obj.svd$center,
                                scale = obj.svd$scale)
# X2 <- sweep(X2, 1, length(keep) / nb_nona[ind2], '*')
proj <- matrix(0, length(ind2), K)
L_all <- U3[ind2, ]  # X2 %*% V
L_all2 <- X2 %*% V
stopifnot(all.equal(L_all, L_all2))
L_norm <- rowSumsSq(L_all)
X2_norm <- rowSumsSq(X2)

for (i in seq_along(ind2)) {
  print(i)

  Q[, dim] <- c(L_all[i, ], sqrt(X2_norm[i] - L_norm[i]))
  R <- crossprod(Q)
  eig <- .Internal(La_rs(R, FALSE))  # eigen_real_symmetric_with_vectors
  PC.aug <- U4 %*% sweep(eig$vectors[, dim - 1:K + 1], 2,
                         sqrt(eig$values[dim - 1:K + 1]), '*')

  proc <- get_procrustes(PC.ref, PC.aug[-nrow(PC.aug), ])
  proj[i, ] <- apply_procrustes(PC.aug[nrow(PC.aug), , drop = FALSE], proc)
}

L <- L_all[i, ]
G <- sqrt(X2_norm[i] - L_norm[i])
d2 <- c(d, 0)
d3 <- d2^2
R2 <- diag(d3)
microbenchmark::microbenchmark(
  EIGEN = {
    Q[, dim] <- c(L, G)
    R <- crossprod(Q)
    eig <- eigen(R, symmetric = TRUE)
    U3 <- U2 %*% eig$vectors[, 1:K]
  },
  EIGEN2 = {
    Q[, dim] <- c(L, G)
    R <- crossprod(Q)
    eig <- .Internal(La_rs(R, FALSE))  # eigen_real_symmetric_with_vectors
    U3 <- U2 %*% eig$vectors[, dim - 1:K + 1]
  },
  EIGEN3 = {
    add <- c(L * d, X2_norm[i])
    R2[, dim] <- add
    R2[dim, ] <- add
    eig <- .Internal(La_rs(R2, FALSE))  # eigen_real_symmetric_with_vectors
    U3 <- U2 %*% eig$vectors[, dim - 1:K + 1]
  },
  RANDOM = {
    LG <- c(L, G)
    svd_aug <- RSpectra::svds(A = function(x, args) x * d2 + x[dim] * LG,
                              Atrans = function(x, args) {
                                res <- x * d2
                                res[dim] <- crossprod(x, LG)
                                res
                              }, k = K, nu = 0, nv = K, dim = c(dim, dim),
                              opts = list(tol = 1e-4))
    U4 <- U2 %*% svd_aug$v
  },
  RANDOM2 = {
    DL <- c(d * L, 0)
    LLGG <- drop(crossprod(L) + G^2)
    svd_aug3 <- RSpectra::eigs_sym(A = function(x, args) {
      res <- x * d3 + x[dim] * DL
      res[dim] <- crossprod(x, DL) + x[dim] * LLGG
      res
    }, k = K, n = dim, opts = list(tol = 1e-4))
    U6 <- U2 %*% svd_aug3$vectors
  },
  SVD = {
    Q[, dim] <- c(L, G)
    svd_aug2 <- svd(Q, nu = 0, nv = K)
    U5 <- U2 %*% svd_aug2$v
  }
)
# Unit: microseconds
#    expr      min        lq      mean   median        uq       max neval
#   EIGEN  673.532  688.5540  704.6623  697.917  720.6475   783.546   100
#  EIGEN2  601.110  618.7680  630.0408  625.215  635.1930   783.201   100
#  EIGEN3  553.379  566.7095  580.8214  575.435  592.9510   629.839   100
#  RANDOM 4435.245 4629.3870 4843.2983 4689.961 4892.0265 11414.592   100
# RANDOM2 2278.564 2341.9825 2661.3136 2382.446 2463.6410 20286.928   100
#     SVD  895.315  918.5145  937.6975  932.878  956.5580  1017.373   100


plot(proj[, 5:10], U3[ind2, 5:10], col = col(proj[, 5:10]))

col <- 5:6 + 24
plot(PC.ref[, col], xlab = paste0("PC", col[1]), ylab = paste0("PC", col[2]))
points(U1[, col], col = "red", pch = 20)
points(proj[, col], col = "blue", pch = 20)
# points(proj2[, col], col = "pink", pch = 20)

k <- 12
plot(U3[ind2, k], proj[, k])
(mylm <- lm(proj[, k] ~ U3[ind2, k] + 0))
abline(mylm, col = "red")
