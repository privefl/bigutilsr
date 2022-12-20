#' Regularization with the graphical lasso
#'
#' Use the graphical lasso algorithm to regularize a square symmetric matrix
#' (e.g. a covariance or correlation matrix) by assuming that its inverse has
#' many zeros.
#'
#' @param mat A square symmetric matrix.
#' @param lambda Strength of regularization. It needs to be scaled with `mat`.
#'   It should also be the maximum difference between the two matrices.
#' @param maxiter_outer Maximum number of iterations of the outer loop.
#'   Default is 200.
#' @param maxiter_lasso Maximum number of iterations of each lasso solver.
#'   Default is 200.
#' @param tol Tolerance for assessing convergence. Default is 1e-4 and it needs
#'   to be scaled with `mat`.
#' @param verbose Whether to print iterations and differences. Default is FALSE.
#'
#' @return The regularized matrix, where the diagonal should be the same and
#'   zeros should be kept as well. It also returns the `lambda` used as an attribute.
#' @export
#'
#' @examples
#' (cov <- cov(iris[1:4]))
#' lambda <- 1 / sqrt(nrow(iris))
#' (cov_regul <- regul_glasso(cov, lambda))
regul_glasso <- function(mat,
                         lambda,
                         maxiter_outer = 200,
                         maxiter_lasso = 200,
                         tol = 1e-4,
                         verbose = FALSE) {
  repeat {

    new_mat <- try(
      glasso(mat, lambda, maxiter_outer, maxiter_lasso, tol, verbose),
      silent = TRUE)

    if (inherits(new_mat, "try-error")) {
      message("Divergence! Increasing lambda (x1.5) and retrying..")
      lambda <- lambda * 1.5
    } else break
  }

  structure(new_mat, lambda = lambda, dimnames = dimnames(mat))
}
