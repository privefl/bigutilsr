context("test-hdpca")

X <- readRDS(system.file("testdata", "three-pops.rds", package = "bigutilsr"))
pop <- rep(1:3, c(143, 167, 207))

N <- 300; M <- ncol(X)
ind <- sample(nrow(X), N)
svd <- svd(X[ind, ])
plot(svd$d^2, log = "xy")

expect_equal(U0 <- sweep(svd$u, 2, svd$d, '*'), X[ind, ] %*% svd$v)
U1 <- X[-ind, ] %*% svd$v
U2 <- pca_adjust(svd$d^2, M, N, test.scores = U1, n.spikes.max = 50)
(shrinkage <- attr(U2, "shrinkage"))

# expect_equal(pca_nspike(svd$d^2, M, N, n.spikes.max = 20)$n.spikes, 3)
expect_gte(length(shrinkage), 3)

# col <- 2:3
# plot(U0[, col])
# points(U1[, col], col = "red", pch = 20)
# points(U2[, col], col = "blue", pch = 20)

ref   <- by(U0[, 1:4], pop[ind],  colMeans)
pred1 <- by(U1[, 1:4], pop[-ind], colMeans)
pred2 <- by(U2[, 1:4], pop[-ind], colMeans)
lapply(seq_along(ref), function(k) {
  expect_lt(crossprod(ref[[k]] - pred2[[k]]), crossprod(ref[[k]] - pred1[[k]]))
})
