#### First example ####
X <- readRDS(system.file("testdata", "three-pops.rds", package = "bigutilsr"))
X <- scale(X)
pop <- rep(1:3, c(143, 167, 207))

N <- 300; M <- ncol(X)
ind <- sample(nrow(X), N)

K <- 30; K2 <- 1 * K
svd <- svd(X[ind, ], nu = K2, nv = K2)
U <- svd$u
d <- head(svd$d, K2)
V <- svd$v

U0 <- svd(rbind(X[ind, ], X[-ind, ][1, ]), nu = K, nv = 0)$u

# `%T*%` <- crossprod
# `%*T%` <- tcrossprod

U2 <- cbind(rbind(U, 0), 0)
U2[nrow(U2), ncol(U2)] <- 1
d2 <- c(d, 0)
d3 <- d2^2
y <- X[-ind, ][1, ]  # should be scaled
L <- crossprod(V, y)
H <- y - V %*% L
H <- H / drop(sqrt(crossprod(H)))
G <- crossprod(y, H)
Q <- cbind(rbind(diag(d), 0), 0)
dim <- nrow(Q)

# library(Matrix)
# Q2 <- Matrix(cbind(rbind(diag(d), 0), 0), sparse = TRUE)

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
#    expr      min        lq      mean   median       uq      max neval
#   EIGEN  719.085  755.5455  851.5371  791.729  852.762 1821.178   100
#  EIGEN2  655.715  689.0475  731.1495  710.551  748.762 1140.547   100
#  RANDOM 3495.966 3590.0990 3953.8574 3744.524 3957.863 9077.300   100
# RANDOM2 2156.622 2284.4350 2603.2374 2383.232 2507.912 8717.175   100
#     SVD 1010.374 1077.3485 1160.1941 1111.605 1184.838 1734.130   100

plot(U4, U3); abline(0, 1, col = "red"); abline(0, -1, col = "red")
plot(U5, U3); abline(0, 1, col = "red"); abline(0, -1, col = "red")
plot(U6, U3); abline(0, 1, col = "red"); abline(0, -1, col = "red")
plot(U3[301, ], U0[301, ]); abline(0, 1, col = "red"); abline(0, -1, col = "red")
cor(abs(U3[301, ]), abs(U0[301, ]))

replicate(1000, {
  Q[, dim] <- c(L, G)
  R <- crossprod(Q)
  eig <- eigen(R, symmetric = TRUE)
  U3 <- U2 %*% eig$vectors[, 1:K]
  NULL
})

replicate(5000, {
  Q[, dim] <- c(L, G)
  R <- crossprod(Q)
  eig <- .Internal(La_rs(R, FALSE))
  U3 <- U2 %*% eig$vectors[, dim - 1:K + 1L]
  NULL
})


