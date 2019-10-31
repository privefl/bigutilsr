################################################################################

F_norm2 <- function(X) {
  sum(rowSumsSq(X))
}

################################################################################

# Y = pXR
get_procrustes <- function(Y, X) {
  X_mean <- colMeans(X)
  X.centered <- sweep(X, 2, X_mean)
  svd <- svd(crossprod(Y, X.centered))
  rho <- sum(svd$d) / F_norm2(X.centered)
  R <- tcrossprod(svd$v, svd$u)
  structure(list(R = R, rho = rho, c = colMeans(Y) - crossprod(rho * X_mean, R)),
            class = "procrustes")
}

################################################################################

#' Predict method
#'
#' Predict method for class `procrustes`.
#'
#' @param object Object of class `procrustes`.
#' @param X New matrix to transform.
#' @param ... Not used.
#'
#' @return `X`, transformed.
#'
#' @export
#' @importFrom stats predict
#'
#' @seealso [procrustes()].
#'
predict.procrustes <- function(object, X, ...) {
  sweep(X %*% (object$rho * object$R), 2, object$c, '+')
}

################################################################################

#' Procrustes transform
#'
#' Procrustes transform Y = pXR (after centering), where p is a scaling
#' coefficient and R is a rotation matrix that minimize ||Y - pXR||_F.
#'
#' @param Y Reference matrix.
#' @param X Matrix to transform (`ncol(X) >= ncol(Y)`).
#' @param n_iter_max Maximum number of iterations. Default is `1000`.
#' @param epsilon_min Convergence criterion. Default is `1e-7`.
#'
#' @return Object of class "procrustes", a list with the following elements:
#'   - `$R`: the rotation matrix to apply to `X`,
#'   - `$rho`: the scaling coefficient to apply to `X`,
#'   - `$c`: the column centering to apply to the resulting matrix,
#'   - `$diff`: the average difference between `Y` and `X` transformed.
#'
#' You can use method `predict()` to apply this transformation to other data.
#'
#' @export
#'
#' @examples
#' A <- matrix(rnorm(200), ncol = 20)
#' B <- matrix(rnorm(length(A)), nrow = nrow(A))
#'
#' proc <- procrustes(B, A)
#' str(proc)
#' plot(B, predict(proc, A)); abline(0, 1, col = "red")
#'
procrustes <- function(Y, X, n_iter_max = 1000, epsilon_min = 1e-7) {

  if (ncol(X) == ncol(Y)) {
    proc <- get_procrustes(Y, X)
    X_new <- predict(proc, X)
  } else if (ncol(X) < ncol(Y)) {
    stop("'X' should have at least the same number of columns as 'Y'.")
  } else {
    Z <- matrix(0, nrow(Y), ncol(X) - ncol(Y))
    for (i in seq_len(n_iter_max)) {
      W <- cbind(Y, Z)
      proc <- get_procrustes(W, X)
      X_new <- predict(proc, X)
      Z_new = X_new[, (ncol(Y) + 1):ncol(X), drop = FALSE]
      eps <- F_norm2(Z - Z_new) / F_norm2(sweep(Z_new, 2, colMeans(Z_new)))
      # print(eps)
      if (eps < epsilon_min) break
      Z <- Z_new
    }
  }

  proc$diff <- F_norm2(X_new[, seq_len(ncol(Y))] - Y) / length(Y)
  proc
}

################################################################################
