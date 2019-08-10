################################################################################

#' Gaussian smoothing
#'
#' @param x Numeric vector.
#' @param size Radius of the smoothing (smaller than half of the length of `x`).
#'
#' @return Numeric vector of the same length as `x`, smoothed.
#' @export
#'
#' @examples
#' (x <- rnorm(10))
#' rollmean(x, 3)
rollmean <- function(x, size) {

  len <- 2 * floor(size) + 1
  if (len >= length(x)) stop("Parameter 'size' is too large.")
  if (size < 0) stop("Parameter 'size' must be positive.")

  lims <- qnorm(range(ppoints(len)))
  weights <- dnorm(seq(lims[1], lims[2], length.out = len))

  roll_mean(x, weights)
}

################################################################################
