context("test-pca-project")

#### First example ####
X <- readRDS(system.file("testdata", "three-pops.rds", package = "bigutilsr"))
pop <- rep(1:3, c(143, 167, 207))

N <- 300; M <- ncol(X)
ind <- sample(nrow(X), N)
svd <- svds(X[ind, ], k = 10)
# plot(svd$d^2, log = "xy")
# hist(svd$d[svd$d < 80]^2, breaks = nclass.scottRob)

expect_equal(U0 <- sweep(svd$u, 2, svd$d, '*'), X[ind, ] %*% svd$v)

X.new <- X[-ind, ]
U1 <- X.new %*% svd$v

proj <- pca_OADP_proj(X.new, loadings = svd$v[, 1:5], sval = svd$d)
expect_equal(proj$simple_proj, U1[, 1:5])
U3 <- proj$OADP_proj

all_maha <- by(U0[, 2:3], pop[ind], function(x) covrob_ogk(as.matrix(x)))
pred1 <- lapply(seq_along(all_maha), function(k) {
  maha <- all_maha[[k]]
  stats::mahalanobis(U1[pop[-ind] == k, 2:3], center = maha$center, cov = maha$cov)
})
pred3 <- lapply(seq_along(all_maha), function(k) {
  maha <- all_maha[[k]]
  stats::mahalanobis(U3[pop[-ind] == k, 2:3], center = maha$center, cov = maha$cov)
})
expect_gt(ks.test(pchisq(unlist(pred3), df = 2, lower.tail = TRUE), "punif")$p.value,
          ks.test(pchisq(unlist(pred1), df = 2, lower.tail = TRUE), "punif")$p.value)


# Augmentation, Decomposition and Procrustes (ADP) transformation
PC.ref <- U0[, 1:3]
X.aug <- rbind(0, X[ind, ])

U4 <- t(sapply(1:nrow(X.new), function(i) {
  # print(i)
  X.aug[1, ] <<- X.new[i, ]
  svd.aug <- svds(X.aug, k = 3, nv = 0)
  PC.aug <- sweep(svd.aug$u, 2, svd.aug$d, '*')
  proc <- procrustes(PC.ref, PC.aug[-1, ])
  predict(proc, PC.aug[1, , drop = FALSE])
}))
expect_equal(U4, U3[, 1:3], tolerance = 1e-2)

U5 <- t(sapply(1:nrow(X.new), function(i) {
  # print(i)
  X.aug[1, ] <<- X.new[i, ]
  svd.aug <- svds(X.aug, k = 6, nv = 0)
  PC.aug <- sweep(svd.aug$u, 2, svd.aug$d, '*')
  proc <- procrustes(PC.ref, PC.aug[-1, ])
  predict(proc, PC.aug[1, , drop = FALSE])[1:3]
}))
expect_equal(U5, U3[, 1:3], tolerance = 1e-3)

# col <- 1:2
# plot(U0[, col])
# points(U1[, col], col = "red", pch = 20)
# points(U2[, col], col = "blue", pch = 20)
# points(U3[, col], col = "green", pch = 20)
# points(U4[, col], col = "orange", pch = 20)
# points(U5[, col], col = "purple", pch = 20)


#### Second example ####
N <- 400; M <- 2000
K <- sample(5:12, 1)
U <- matrix(0, N, K); U[] <- rnorm(length(U))
V <- matrix(0, M, K); V[] <- rnorm(length(V))
# X = U V^T + E
X <- tcrossprod(U, V) + 15 * rnorm(N * M)
svd <- svd(scale(X))
# plot(head(svd$d^2, -1), log = "xy", main = K)
expect_lte(abs(pca_nspike(svd$d^2) - K), 2)
