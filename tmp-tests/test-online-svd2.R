ADD <- c(0, 2, 5, 10, 20, 40, 80)
corr <- sapply(ADD, function(add) {

  K <- 30; K2 <- K + add
  svd <- svd(X[ind, ], nu = K2, nv = K2)
  U <- svd$u
  d <- head(svd$d, K2)
  V <- svd$v

  U0 <- svd(rbind(X[ind, ], X[-ind, ][1, ]), nu = K, nv = 0)$u

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

  Q[, dim] <- c(L, G)
  R <- crossprod(Q)
  eig <- .Internal(La_rs(R, FALSE))  # eigen_real_symmetric_with_vectors
  U3 <- U2 %*% eig$vectors[, dim - 1:K + 1]

  cor(abs(U3[301, ]), abs(U0[301, ]))
})
plot(ADD, corr, ylim = c(min(corr), 1))
