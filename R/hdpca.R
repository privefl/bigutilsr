################################################################################

#' Number of spikes in PCA
#'
#' Estimate the number of distant spikes based on the histogram of eigenvalues.
#'
#' @param eigval Eigenvalues (squared singular values).
#' @param nboot Number of bootstrap replicates to estimate limits more robustly.
#'   Default is `100`.
#' @inheritParams hist_out
#'
#' @return The estimated number of distant spikes.
#'
#' @seealso [hdpca::select.nspike()]
#'
#' @export
#'
#' @examples
#' N <- 400; M <- 2000; K <- 8
#' U <- matrix(0, N, K); U[] <- rnorm(length(U))
#' V <- matrix(0, M, K); V[] <- rnorm(length(V))
#' # X = U V^T + E
#' X <- tcrossprod(U, V) + 15 * rnorm(N * M)
#' pca <- prcomp(X)
#' eigval <- pca$sdev^2
#' plot(head(eigval, -1), log = "xy", pch = 20)
#' pca_nspike(eigval)
#'
pca_nspike <- function(eigval, breaks = "FD", nboot = 100) {
  lim_up <- hist_out(eigval, breaks = "FD", nboot = nboot)$lim[2]
  sum(eigval > lim_up)
}

################################################################################

#' @inherit hdpca::pc_adjust title description details params references
#'
#' @return
#' A matrix containing the bias-adjusted PC scores.
#' The dimension of the matrix is the same as the dimension of `test.scores`.
#'
#' Also, an attribute `attr(*, "shrinkage")` containing the shrinkage factors.
#' Note that the number of shrinkage factors can be smaller than the number of
#' columns of `test.scores`; it corresponds to the estimated number of spikes.
#'
#' @seealso [hdpca::pc_adjust()]
#'
#' @export
#'
#' @examples
#' X <- readRDS(system.file("testdata", "three-pops.rds", package = "bigutilsr"))
#' N <- 400; M <- ncol(X)
#' ind <- sample(nrow(X), N)
#' # Compute SVD using one part of samples
#' svd <- svd(X[ind, ])
#' U <- sweep(svd$u, 2, svd$d, '*')
#' plot(U[, 2:3])
#' # Projecting other samples
#' U1 <- X[-ind, ] %*% svd$v
#' points(U1[, 2:3], col = "red", pch = 20)               # shrunk towards 0
#' U2 <- pca_adjust(U1, svd$d^2, M, N)
#' points(U2[, 2:3], col = "blue", pch = 20)              # unshrunk
#' attr(U2, "shrinkage")
#'
pca_adjust <- function(test.scores, train.eval, p,
                       n = length(train.eval),
                       method = c("d.gsp", "l.gsp", "osp"),
                       n.spikes = pca_nspike(train.eval)) {

  stopifnot(length(dim(test.scores)) == 2)

  shrinkage <- hdpca::hdpc_est(train.eval, p, n, method,
                               n.spikes = n.spikes)$shrinkage

  m <- min(length(shrinkage), ncol(test.scores))
  for (k in seq_len(m)) test.scores[, k] <- test.scores[, k] / shrinkage[k]

  structure(test.scores, shrinkage = shrinkage[1:m])
}

################################################################################
