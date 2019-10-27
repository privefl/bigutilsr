################################################################################

#' OADP projection
#'
#' Online Augmentation, Decomposition, and Procrustes (OADP) projection of
#' PC loadings onto some study data `X`.
#'
#' @param X Data to get PC loadings into.
#' @param loadings PC loadings of the reference PCA to project.
#' @param sval Singular values of the reference PCA (sqrt of the eigen values).
#'   Only the `ncol(loadings)` first ones will be used.
#'
#' @return
#'  - `pca_OADP_proj()`: A list with the simple projection `X %*% loadings`
#'    and the projection based on OADP.
#' @export
#'
#' @examples
#' X <- readRDS(system.file("testdata", "three-pops.rds", package = "bigutilsr"))
#' N <- 400; M <- ncol(X)
#' ind <- sample(nrow(X), N)
#' # Compute SVD using one part of samples
#' svd <- svd(X[ind, ])
#' U <- sweep(svd$u, 2, svd$d, '*')
#' col <- 2:3
#' plot(U[, col])
#' points(cbind(0, 0), pch = 8, col = "green", cex = 2)
#' # Projecting other samples
#' proj <- pca_OADP_proj(X = X[-ind, ], loadings = svd$v[, 1:5], sval = svd$d)
#' points(proj$simple_proj[, col], col = "red", pch = 20)     # shrunk towards 0
#' points(proj$OADP_proj[, col], col = "blue", pch = 20)      # unshrunk
#'
pca_OADP_proj <- function(X, loadings, sval) {

  XV     <- X %*% loadings
  X_norm <- rowSumsSq(X)

  list(simple_proj = XV, OADP_proj = pca_OADP_proj2(XV, X_norm, sval))
}

################################################################################

#' @rdname pca_OADP_proj
#'
#' @param XV `X %*% loadings`
#' @param X_norm Vector of sums of squared rows (e.g. `rowSums(X^2)`).
#'
#' @return
#'  - `pca_OADP_proj2()`: The projection based on OADP only
#'    (a matrix of same size of `XV`).
#' @export
#'
pca_OADP_proj2 <- function(XV, X_norm, sval) {

  m <- nrow(XV)
  K <- ncol(XV)
  d <- sval[1:K]

  QtQ <- diag(c(d^2, 0))

  proj <- matrix(0, m, K)

  for (i in seq_len(m)) {

    QtQ[K + 1, ] <- QtQ[, K + 1] <- c(XV[i, ] * d, X_norm[i])
    eig <- eigen(QtQ, symmetric = TRUE)

    V2 <- sweep(eig$vectors, 2, sqrt(eig$values), '*')[, 1:K, drop = FALSE]
    svd2 <- svd(sweep(V2[1:K, , drop = FALSE], 1, d, '*'))
    rho <- sum(svd2$d) / cumsum(rowSumsSq(V2))[K]
    proj[i, ] <- rho * tcrossprod(
      V2[K + 1, , drop = FALSE] %*% svd2$v, svd2$u)
  }

  proj
}

################################################################################
