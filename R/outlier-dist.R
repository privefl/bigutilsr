################################################################################

#' Transform matrix
#'
#' Transform matrix to use Mahalanobis distance instead of Euclidean one.
#'
#' @param U A matrix (e.g. PC scores).
#' @param estim List of location and scatter estimates, `$cov` and `$center`.
#'
#' @return `U`, transformed.
#' @export
#'
#' @examples
#' X <- readRDS(system.file("testdata", "three-pops.rds", package = "bigutilsr"))
#' svd <- svds(scale(X), k = 5)
#'
#' U <- svd$u
#' dist1 <- dist_ogk(U)
#'
#' U.maha <- maha_trans(U)
#' dist2 <- rowSums(U.maha^2)
#' all.equal(dist2, dist1)
#'
maha_trans <- function(U, estim = covrob_ogk(U)) {

  eigs <- eigen(unname(estim$cov), symmetric = TRUE)
  sqrt_inv_cov <- sweep(eigs$vectors, 2, sqrt(eigs$values), '/')

  sweep(U, 2, estim$center, '-') %*% sqrt_inv_cov
}

################################################################################

#' Local Outlier Factor (LOF)
#'
#' LOF: Identifying Density-Based Local Outliers.
#'
#' @param U A matrix, from which to detect outliers (rows). E.g. PC scores.
#' @param seq_k Sequence of numbers of nearest neighbors to use.
#'   If multiple `k` are provided, this returns the combination of statistics.
#'   Default is `c(4, 10, 30)` and use `max` to combine (see `combine`).
#' @param combine How to combine results for multiple `k`? Default uses `max`.
#' @param robMaha Whether to use a robust Mahalanobis distance instead of the
#'   normal euclidean distance? Default is `FALSE`, meaning using euclidean.
#' @param log Whether to return the logarithm of LOFs? Default is `TRUE`.
#' @param ncores Number of cores to use. Default is `1`.
#'
#' @seealso [prob_dist()]
#'
#' @references
#' Breunig, Markus M., et al. "LOF: identifying density-based local outliers."
#' ACM sigmod record. Vol. 29. No. 2. ACM, 2000.
#'
#' @export
#'
#' @examples
#' X <- readRDS(system.file("testdata", "three-pops.rds", package = "bigutilsr"))
#' svd <- svds(scale(X), k = 10)
#'
#' llof <- LOF(svd$u)
#' hist(llof, breaks = nclass.scottRob)
#' tukey_mc_up(llof)
#'
#' llof_maha <- LOF(svd$u, robMaha = TRUE)
#' hist(llof_maha, breaks = nclass.scottRob)
#' tukey_mc_up(llof_maha)
#'
#' lof <- LOF(svd$u, log = FALSE)
#' hist(lof, breaks = nclass.scottRob)
#' str(hist_out(lof))
#' str(hist_out(lof, nboot = 100))
#' str(hist_out(lof, nboot = 100, breaks = "FD"))
#'
LOF <- function(U, seq_k = c(4, 10, 30), combine = max,
                robMaha = FALSE, log = TRUE, ncores = 1) {

  if (robMaha) U <- maha_trans(U)

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
#' @param kNN Number of nearest neighbors to use. Default is `5`.
#'
#' @seealso [LOF()]
#'
#' @references
#' Kriegel, Hans-Peter, et al. "LoOP: local outlier probabilities." Proceedings
#' of the 18th ACM conference on Information and knowledge management. ACM, 2009.
#'
#' @export
#'
#' @examples
#' X <- readRDS(system.file("testdata", "three-pops.rds", package = "bigutilsr"))
#' svd <- svds(scale(X), k = 10)
#' U <- svd$u
#'
#' test <- prob_dist(U)
#' plof <- test$dist.self / test$dist.nn
#' plof_ish <- test$dist.self / sqrt(test$dist.nn)
#' plot(U[, 1:2], col = (plof_ish > tukey_mc_up(plof_ish)) + 1, pch = 20)
#' plot(U[, 3:4], col = (plof_ish > tukey_mc_up(plof_ish)) + 1, pch = 20)
#' plot(U[, 5:6], col = (plof_ish > tukey_mc_up(plof_ish)) + 1, pch = 20)
#'
prob_dist <- function(U, kNN = 5, robMaha = FALSE, ncores = 1) {

  if (robMaha) U <- maha_trans(U)

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

#' Geometric median
#'
#' Compute the geometric median, i.e. the point that minimizes the sum of all
#' Euclidean distances to the observations (rows of `U`).
#'
#' @param U A matrix (e.g. PC scores).
#' @param tol Convergence criterion. Default is `1e-10`.
#' @param maxiter Maximum number of iterations. Default is `1000`.
#' @param by_grp Possibly a vector for splitting rows of `U` into groups before
#'   computing the geometric mean for each group. Default is `NULL` (ignored).
#'
#' @return The geometric median of all rows of `U`, a vector of the same size
#'   as `ncol(U)`. If providing `by_grp`, then a matrix with rows being the
#'   geometric median within each group.
#' @export
#'
#' @examples
#' X <- readRDS(system.file("testdata", "three-pops.rds", package = "bigutilsr"))
#' pop <- rep(1:3, c(143, 167, 207))
#'
#' svd <- svds(scale(X), k = 5)
#' U <- sweep(svd$u, 2, svd$d, '*')
#' plot(U, col = pop, pch = 20)
#'
#' med_all <- geometric_median(U)
#' points(t(med_all), pch = 20, col = "blue", cex = 4)
#'
#' med_pop <- geometric_median(U, by_grp = pop)
#' points(med_pop, pch = 20, col = "blue", cex = 2)
#'
geometric_median <- function(U, tol = 1e-10, maxiter = 1000, by_grp = NULL) {

  if (!is.null(by_grp))
    return(do.call("rbind", by(U, by_grp, geometric_median)))

  # Weiszfeld's algorithm
  u.old <- colMeans(U)

  for (k in seq_len(maxiter)) {
    norm <- sqrt(rowSums(sweep(U, 2, u.old, '-')^2))
    u.new <- colSums(sweep(U, 1, norm, '/')) / sum(1 / norm)
    diff <- max(abs(u.new - u.old))
    if (diff < tol) break
    u.old <- u.new
  }

  if (k == maxiter)
    warning("The maximum number of iterations has been reached.")

  u.new
}

################################################################################
