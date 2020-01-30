################################################################################

get_sigma0 <- function(x, tmp_dev, tmp_med)
  scaleTau2_vector_rcpp(x, tmp_dev, tmp_med)[2]

ogk_step_r <- function(X) {

  p <- ncol(X)
  S <- matrix(1, p, p)
  X.col.scaled <- vector("list", p)

  tmp_col <- double(nrow(X))
  tmp_dev <- double(nrow(X))
  tmp_med <- double(nrow(X))

  for (j in seq_len(p)) {
    X_j <- X[, j]
    sigma0_j <- get_sigma0(X_j, tmp_dev, tmp_med)
    col_j <- X.col.scaled[[j]] <- X_j / sigma0_j
    for (k in seq_len(j - 1)) {
      col_k <- X.col.scaled[[k]]
      sigma0_sum  <- get_sigma0(sum_in_temp(col_j, col_k, tmp_col), tmp_dev, tmp_med)
      sigma0_diff <- get_sigma0(sub_in_temp(col_j, col_k, tmp_col), tmp_dev, tmp_med)
      S[j, k] <- S[k, j] <- (sigma0_sum ** 2 - sigma0_diff ** 2) / 4
    }
  }

  if (anyNA(S)) {
    bad_col <- unique(which(is.na(S) & upper.tri(S), arr.ind = TRUE)[, 2])
    stop2("Problem with columns [%s] of 'U'.", paste(bad_col, collapse = ", "))
  }

  do.call("cbind", X.col.scaled) %*% eigen(S, symmetric = TRUE)$vectors
}

################################################################################

#' @inherit rrcov::CovOgk title description details references
#'
#' @param U A matrix with no missing values and at least 2 columns.
#' @param niter Number of number of iterations for the first step of the algorithm,
#'   usually 1 or 2 since iterations beyond the second do not lead to improvement.
#' @param beta Coverage parameter for the final reweighted estimate.
#'   Default is `0.9`.
#'
#' @return `covrob_ogk()`: list of robust estimates, `$cov` and `$center`.
#' @export
#'
#' @seealso [rrcov::CovOgk()]
#'
covrob_ogk <- function(U, niter = 2, beta = 0.9) {

  bigassertr::assert_class(U, "matrix")
  bigassertr::assert_nona(U)

  Z <- U
  for (k in seq_len(niter)) {
    Z <- ogk_step_r(Z)
  }

  d <- dist_scaleTau2_matrix_rcpp(Z)

  df <- ncol(U)
  cdelta <- median(d) / qchisq(0.5, df)
  quantile <- qchisq(beta, df)

  U.in <- U[d < (quantile * cdelta), ]
  wcenter <- colMeans(U.in)
  wcov <- crossprod(sweep(U.in, 2, wcenter, '-')) / nrow(U.in)

  list(cov = wcov, center = wcenter)
}

################################################################################

#' @rdname covrob_ogk
#'
#' @return `dist_ogk()`: vector of robust Mahalanobis (squared) distances.
#' @export
#'
#' @seealso [stats::mahalanobis()]
#'
#' @examples
#' X <- readRDS(system.file("testdata", "three-pops.rds", package = "bigutilsr"))
#' svd <- svds(scale(X), k = 5)
#'
#' U <- svd$u
#' dist <- dist_ogk(U)
#' str(dist)
dist_ogk <- function(U, niter = 2, beta = 0.9) {
  estim <- covrob_ogk(U, niter, beta)
  stats::mahalanobis(U, estim$center, estim$cov)
}

################################################################################
