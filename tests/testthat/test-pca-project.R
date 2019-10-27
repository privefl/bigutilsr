context("test-pca-project")

#### First example ####
X <- readRDS(system.file("testdata", "three-pops.rds", package = "bigutilsr"))
pop <- rep(1:3, c(143, 167, 207))

N <- 300; M <- ncol(X)
ind <- sample(nrow(X), N)
svd <- svd(X[ind, ])
# plot(svd$d^2, log = "xy")
# hist(svd$d[svd$d < 80]^2, breaks = nclass.scottRob)

expect_equal(U0 <- sweep(svd$u, 2, svd$d, '*'), X[ind, ] %*% svd$v)
U1 <- X[-ind, ] %*% svd$v

proj <- pca_OADP_proj(X[-ind, ], loadings = svd$v[, 1:5], sval = svd$d)
expect_equal(proj$simple_proj, U1[, 1:5])
U3 <- proj$OADP_proj

# col <- 2:3
# plot(U0[, col])
# points(U1[, col], col = "red", pch = 20)
# points(U2[, col], col = "blue", pch = 20)
# points(U3[, col], col = "green", pch = 20)

ref   <- by(U0[, 1:3], pop[ind],  colMeans)
pred1 <- by(U1[, 1:3], pop[-ind], colMeans)
pred3 <- by(U3[, 1:3], pop[-ind], colMeans)
lapply(seq_along(ref), function(k) {
  expect_lt(crossprod(ref[[k]] - pred3[[k]]), crossprod(ref[[k]] - pred1[[k]]))
})


#### Second example ####
N <- 400; M <- 2000
K <- sample(5:12, 1)
U <- matrix(0, N, K); U[] <- rnorm(length(U))
V <- matrix(0, M, K); V[] <- rnorm(length(V))
# X = U V^T + E
X <- tcrossprod(U, V) + 15 * rnorm(N * M)
svd <- svd(scale(X))
# plot(head(svd$d^2, -1), log = "xy", main = K)
expect_equal(pca_nspike(svd$d^2), K)
