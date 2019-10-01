library(bigsnpr)
bedfile <- download_1000G("tmp-data")
(obj.bed <- bed(bedfile))

N <- 1800; M <- ncol(obj.bed)
set.seed(1); ind <- sample(nrow(obj.bed), N)
K <- 30
obj.svd <- bed_autoSVD2(obj.bed, ind.row = ind, k = K, ncores = nb_cores())
U <- obj.svd$u
d <- obj.svd$d
V <- obj.svd$v

keep <- attr(obj.svd, "subset.col")
{
  U2 <- big_apply(obj.bed, function(X, ind, ind.col, center, scale, V) {
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

PC.ref <- predict(obj.svd)


# U4 <- cbind(rbind(U, 0), 0)
# U4[nrow(U4), ncol(U4)] <- 1
R <- diag(c(d^2, 0))

set.seed(2); ind2 <- rows_along(obj.bed)[-ind] #sample(rows_along(obj.bed)[-ind], 200)
X2 <- bigsnpr:::read_bed_scaled(obj.bed, ind2, keep,
                                center = obj.svd$center,
                                scale = obj.svd$scale)
# X2 <- sweep(X2, 1, length(keep) / nb_nona[ind2], '*')
proj <- matrix(0, length(ind2), K)
L_all <- U3[ind2, ]  # X2 %*% V
L_all2 <- X2 %*% V
stopifnot(all.equal(L_all, L_all2))
Rcpp::sourceCpp('tmp-tests/rowSumsSq.cpp')
# L_norm <- rowSumsSq(L_all)
X2_norm <- rowSumsSq(X2)

for (i in seq_along(ind2)) {

  print(i)

  R[K + 1, ] <- R[, K + 1] <- c(L_all[i, ] * d, X2_norm[i])

  eig <- .Internal(La_rs(R, FALSE))  # eigen_real_symmetric_with_vectors
  # PC.aug <- U4 %*% sweep(eig$vectors, 2, sqrt(eig$values), '*')[, dim - 1:K + 1]
  # d2 <- eig$values[dim - 1:K + 1]
  # VD <- sweep(, 2, sqrt(d2), '*')
  # PC.aug2 <- U %*% eig$vectors[1:K, (K + 1):2]
  # stopifnot(all.equal(PC.aug2, PC.aug[-nrow(PC.aug), ]))

  # proc <- get_procrustes(PC.ref, PC.aug[-nrow(PC.aug), ])
  # proj[i, ] <- apply_procrustes(PC.aug[nrow(PC.aug), , drop = FALSE], proc)

  # X <- PC.aug2
  # X_mean <- colMeans(X)
  # X.centered <- PC.aug2 #sweep(X, 2, X_mean)
  # Y <- PC.ref
  # X.norm = sum(X.centered ** 2)
  # all.equal(diag(crossprod(PC.aug2)), d2)
  # diag(crossprod(PC.aug2)) / d2
  # X.norm <- sum(d2)
  # stopifnot(all.equal(crossprod(Y, X.centered), sweep(VD[1:K, ], 1, d, '*')))
  svd <- svd(eig$vectors[1:K, (K + 1):2])
  rho <- sum(svd$d) / drop(crossprod(svd$d))
  proj[i, ] <- rho * tcrossprod(
    eig$vectors[K + 1, (K + 1):2, drop = FALSE] %*% svd$v, svd$u)
}

# microbenchmark::microbenchmark(
#   sum(rowSumsSq(PC.aug2)),
#   norm(PC.aug2, type = "F")^2
# )
proj2 <- sweep(proj, 2, d, '*')

plot(proj[, 5:10], U3[ind2, 5:10], col = col(proj[, 5:10]))

col <- 5:6 + 10
plot(PC.ref[, col], xlab = paste0("PC", col[1]), ylab = paste0("PC", col[2]))
points(U3[ind2, col], col = "red", pch = 20)
points(proj2[, col], col = "blue", pch = 20)
