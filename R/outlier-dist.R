################################################################################

#' @importFrom robust covRob
#' @export
robust::covRob

################################################################################

#' Local Outlier Factor (LOF)
#'
#' LOF: Identifying Density-Based Local Outliers.
#'
#' @param U A matrix, from which to detect outliers (rows). E.g. PC scores.
#' @param seq_k Sequence of numbers of nearest neighbours to use.
#'   If multiple `k` are provided, this returns the combination of statistics.
#'   Default is `c(4, 10, 30)` and use `max` to combine (see `combine`).
#' @param combine How to combine results for multiple `k`? Default uses `max`.
#' @param robMaha Whether to use a robust Mahalanobis distance instead of the
#'   normal euclidean distance? Default is `FALSE`, meaning using euclidean.
#' @param log Whether to return the logarithm of LOFs? Default is `TRUE`.
#' @param ncores Number of cores to use. Default is `1`.
#'
#' @references
#' Breunig, Markus M., et al. "LOF: identifying density-based local outliers."
#' ACM sigmod record. Vol. 29. No. 2. ACM, 2000.
#'
#' @export
#'
#' @examples
#' X <- readRDS(system.file("testdata", "three-pops.rds", package = "bigutilsr"))
#' svd <- svd(scale(X))
#'
#' llof <- LOF(svd$u[, 1:10])
#' hist(llof, breaks = nclass.scottRob)
#' tukey_mc_up(llof)
#'
#' lof <- LOF(svd$u[, 1:10], log = FALSE)
#' hist(lof, breaks = nclass.scottRob)
#' str(hist_out(lof))
#' str(hist_out(lof, nboot = 100))
#' str(hist_out(lof, nboot = 100, breaks = "FD"))
#'
LOF <- function(U, seq_k = c(4, 10, 30), combine = max,
                robMaha = FALSE, log = TRUE, ncores = 1) {

  if (robMaha) {
    maha <- covRob(U, estim = "pairwiseGK", distance = FALSE, corr = FALSE)
    eigs <- eigen(unname(maha$cov), symmetric = TRUE)
    U <- U %*% sweep(eigs$vectors, 2, sqrt(eigs$values), '/')
  }

  knn <- knn_parallel(U, k = max(seq_k) + 1, ncores = ncores)
  ids <- knn$nn.idx[, -1, drop = FALSE]
  dists <- knn$nn.dists[, -1, drop = FALSE]

  lof3 <- sapply(seq_k, function(k) {

    lrd <- sapply(seq_len(nrow(dists)), function(i) {
      maxs <- pmax(dists[ids[i, 1:k], k], dists[i, 1:k])
      1 / mean(maxs)
    })

    lof.num <- sapply(seq_len(nrow(ids)), function(i) {
      mean(lrd[ids[i, 1:k]])
    })

    lof.num / lrd
  })

  lof.comb <- apply(lof3, 1, combine)

  `if`(log, log(lof.comb), lof.comb)
}

################################################################################

#' Probabilistic set distance
#'
#' @inheritParams LOF
#' @param kNN Number of nearest neighbours to use. Default is `5`.
#'
#' @references
#' Kriegel, Hans-Peter, et al. "LoOP: local outlier probabilities." Proceedings
#' of the 18th ACM conference on Information and knowledge management. ACM, 2009.
#'
#' @export
#'
#' @examples
#' X <- readRDS(system.file("testdata", "three-pops.rds", package = "bigutilsr"))
#' svd <- svd(scale(X))
#' U <- svd$u[, 1:10]
#'
#' test <- prob_dist(U)
#' plof <- test$dist.self / test$dist.nn
#' plof_ish <- test$dist.self / sqrt(test$dist.nn)
#' plot(U[, 1:2], col = (plof_ish > tukey_mc_up(plof_ish)) + 1, pch = 20)
#' plot(U[, 3:4], col = (plof_ish > tukey_mc_up(plof_ish)) + 1, pch = 20)
#' plot(U[, 5:6], col = (plof_ish > tukey_mc_up(plof_ish)) + 1, pch = 20)
#'
prob_dist <- function(U, kNN = 5, robMaha = FALSE, ncores = 1) {

  if (robMaha) {
    maha <- covRob(U, estim = "pairwiseGK", distance = FALSE, corr = FALSE)
    eigs <- eigen(unname(maha$cov), symmetric = TRUE)
    U <- U %*% sweep(eigs$vectors, 2, sqrt(eigs$values), '/')
  }

  knn <- knn_parallel(U, k = kNN + 1, ncores = ncores)
  ids <- knn$nn.idx[, -1, drop = FALSE]
  dists <- knn$nn.dists[, -1, drop = FALSE]

  plof.num <- sqrt(rowMeans(dists^2))
  plof.deno <- sapply(bigparallelr::rows_along(ids), function(i) {
    mean(plof.num[ids[i, ]])
  })

  list(dist.self = plof.num, dist.nn = plof.deno)
}

################################################################################