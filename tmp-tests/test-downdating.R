set.seed(1)
n <- 2e3
p <- 20e3
U <- matrix(0, n, 10); U[] <- rnorm(length(U))
V <- matrix(0, p, 10); V[] <- rnorm(length(V))
X <- tcrossprod(U, V) + 5 * matrix(rnorm(n * p), n, p)
system.time(svd2 <- PRIMME::svds(X, 10, isreal = TRUE))

ind <- -(1:100)
A <- -X[, -ind]
UT_A <- crossprod(svd2$u, A)
M_A <- A - svd2$u %*% UT_A
svd_A <- svd(M_A)
BT_V <- svd2$v[-ind, ]
# B <- matrix(0, ncol(X), length(ind)); B[cbind(-ind, seq_along(ind))] <- 1
# all.equal(BT_V, crossprod(B, svd2$v))
M_B <- -tcrossprod(svd2$v, BT_V)
M_B[cbind(-ind, seq_along(ind))] <- M_B[cbind(-ind, seq_along(ind))] + 1
svd_B <- svd(M_B)

Q1 <- UT_A %*% BT_V + diag(svd2$d)
VD_A <- sweep(svd_A$v, 2, svd_A$d, '*')
Q2 <- crossprod(VD_A, BT_V)
VD_B <- sweep(svd_B$v, 2, svd_B$d, '*')
Q3 <- UT_A %*% VD_B
Q4 <- crossprod(VD_A, VD_B)
Q <- cbind(rbind(Q1, Q2), rbind(Q3, Q4))
svd_Q <- svd(Q, nu = 10, nv = 10)
U2 <- cbind(svd2$u, svd_A$u) %*% svd_Q$u
V2 <- cbind(svd2$v, svd_B$u) %*% svd_Q$v

# MTM <- crossprod(A) - crossprod(crossprod(svd2$u, A))

X2 <- X[, ind]
u0 <- svd2$u
system.time(
  svd2.1 <- PRIMME::svds(X2, 10, isreal = TRUE))
system.time(
  svd2.2 <- PRIMME::svds(X2, 10, isreal = TRUE, u0 = u0, v0 = V[ind, ]))
system.time(
  svd2.3 <- PRIMME::svds(X2, 10, isreal = TRUE, u0 = u0, v0 = svd(V[ind, ])$u))
system.time(
  svd2.4 <- PRIMME::svds(X2, 10, isreal = TRUE, u0 = u0,
                         v0 = crossprod(X2, sweep(u0, 2, svd2$d, '/'))))
system.time(
  svd2.5 <- PRIMME::svds(X2, 10, isreal = TRUE, u0 = NULL, v0 = V2[ind, ]))
system.time(
  svd2.6 <- PRIMME::svds(X2, 10, isreal = TRUE, u0 = U2, v0 = NULL))

plot(U2, svd2.1$u); abline(0, 1, col = "red"); abline(0, -1, col = "red")
plot(V2[ind, ], svd2.1$v); abline(0, 1, col = "red"); abline(0, -1, col = "red")
all.equal(svd_Q$d[1:10], svd2.1$d)
