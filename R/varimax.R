################################################################################

#' Varimax rotation
#'
#' @param X A matrix with more rows than columns.
#' @param normalize Whether to apply Kaiser normalization? See [stats::varimax].
#'   Default is `FALSE`.
#' @param reorder Whether to permute rotation vectors to maximize the conservation
#'   of the order of the initial columns of `X`. Default is `TRUE`.
#' @param rotmat Whether to return the rotation matrix `rot`, or the rotated
#'   matrix `X %*% rot` (the default, `FALSE`).
#'
#' @return Either the rotation matrix `rot`, or the rotated matrix `X %*% rot`,
#'   depending on `rotmat`.
#' @export
#'
#' @examples
#' X <- as.matrix(iris[1:4])
#' X_rot <- varimax2(X)
#' X_rot2 <- varimax(X, normalize = FALSE)$loadings[]
#' all.equal(X_rot2, X_rot[, c(3, 2, 1, 4)], check.attributes = FALSE)
#' varimax2(X, rotmat = TRUE)
#'
#' X2 <- prcomp(X)$x
#' X2_rot <- varimax2(X2)
#' X2_rot2 <- varimax(X2, normalize = FALSE)$loadings[]
#' all.equal(X2_rot, X2_rot2, check.attributes = FALSE)
varimax2 <- function(X, normalize = FALSE, reorder = TRUE, rotmat = FALSE) {

  rot <- stats::varimax(X, normalize = normalize)$rotmat

  if (reorder) {

    ord <- integer(ncol(X))
    ind_keep <- 1:ncol(X)
    for (k in 1:ncol(X)) {
      ind <- which.max(colSums(rot[1:k, ind_keep, drop = FALSE]^2))
      ord[k]   <- ind_keep[ind]
      ind_keep <- ind_keep[-ind]
    }

    rot <- rot[, ord]
  }

  if (rotmat) rot else X %*% rot
}

################################################################################
