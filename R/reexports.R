################################################################################

#' Deprecated
#'
#' @param data A matrix.
#' @param ... Not used.
#'
#' @seealso [bigutilsr::covrob_ogk()] [bigutilsr::dist_ogk()]
#' @export
covRob <- function(data, ...) {
  .Deprecated("bigutilsr::covrob_ogk()' or 'bigutilsr::dist_ogk()")
  res <- covrob_ogk(data)
  res$dist <- stats::mahalanobis(data, res$center, res$cov)
  res
}

################################################################################

#' @importFrom RSpectra svds
#' @export
RSpectra::svds

################################################################################
