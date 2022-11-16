################################################################################

#' Outlier detection threshold (upper)
#'
#' Outlier detection threshold (upper) based on Tukey's rule, corrected for
#' skewness using the 'medcouple', and possibly corrected for multiple testing.
#'
#' @param x Numeric vector. Should be somewhat normally distributed.
#' @inheritParams robustbase::adjboxStats
#' @param coef number determining how far 'whiskers' extend out from the box.
#'   If `NULL` (default), this is computed to get an type-I error of `alpha`,
#'   after adjusting for multiple testing. A standard value to use is `1.5`.
#' @param alpha See `coef`. Default is `0.05`.
#'
#' @inherit robustbase::adjbox references
#'
#' @seealso [robustbase::adjbox()]
#' @export
#'
#' @examples
#' hist(x <- c(rnorm(3, m = 6), rnorm(1e4, m = 0)))
#' (q <- tukey_mc_up(x))
#' abline(v = q, col = "red")
#' which(x > q)
tukey_mc_up <- function(x, coef = NULL, alpha = 0.05, a = -4, b = 3) {

  q <- Inf

  repeat {

    x <- x[x < q]

    if (is.null(coef)) {
      m <- sum(!is.na(x))
      coef <- (qnorm(log(1 - alpha) / m, log.p = TRUE) - qnorm(0.75)) /
        (qnorm(0.75) - qnorm(0.25))
    }

    q.new <- robustbase::adjboxStats(
      x, coef = coef, a = a, b = b, doReflect = FALSE, doScale = FALSE,
      do.conf = FALSE, do.out = FALSE)$fence[2]

    if (q.new == q) break
    q <- q.new
  }

  q
}

################################################################################
