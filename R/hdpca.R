################################################################################

#' Number of spikes in PCA
#'
#' Estimate the number of distant spikes based on the histogram of eigenvalues.
#'
#' @param eigval Eigenvalues (squared singular values).
#' @param nb_in_bin Mean number of values in each bin of the histogram.
#'
#' @return The estimated number of distant spikes.
#'
#' @seealso [hdpca::select.nspike()]
#'
#' @export
#'
pca_nspike <- function(eigval, nb_in_bin = 10) {

  S <- utils::head(sort(eigval, decreasing = TRUE), -1)

  q <- Inf
  repeat {
    # histogram without outliers
    ind_no_out <- which(S < q)
    nbreak <- round(length(ind_no_out) / nb_in_bin)
    h <- graphics::hist(S[ind_no_out], breaks = nbreak, plot = FALSE)
    # find empty regions on the right
    L <- length(h$mids)
    wh0 <- which(h$counts == 0)
    ind_mid <- wh0[wh0 > L/3][1]
    if (is.na(ind_mid)) break
    # define new threshold to refine histogram
    q <- h$mids[ind_mid]
  }

  sum(S > q)
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
