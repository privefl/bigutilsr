################################################################################

#' @inherit grDevices::nclass.scott title params return
#'
#' @export
#'
#' @references
#' Scott, D. W. (1979). On optimal and data-based histograms.
#' Biometrika, 66, 605–610. doi: 10.2307/2335182.
#'
#' @examples
#' x <- rnorm(1000)
#' hist(x, breaks = nclass.scott)
#' hist(x, breaks = nclass.scottRob)
#'
#' x2 <- c(x, rnorm(50, mean = 50))
#' hist(x2, breaks = nclass.scott)
#' hist(x2, breaks = nclass.scott,    xlim = c(-5, 5))
#' hist(x2, breaks = nclass.scottRob, xlim = c(-5, 5))
#'
nclass.scottRob <- function(x) {
  h <- 3.5 * mad(x) * length(x)^(-1/3)
  if (h > 0) max(1, ceiling(diff(range(x)) / h)) else 1L
}

################################################################################

#' Outlier detection (histogram)
#'
#' Outlier detection based on departure from histogram.
#' Suitable for compact values (need a space between main values and outliers).
#'
#' @param x Numeric vector (with compact values).
#' @param breaks Same parameter as for `hist()`. Default uses a robust version
#'   of Scott's rule. You can also use `"FD"` or `nclass.FD` for a bit more bins.
#' @param pmax_out Percentage at each side that can be considered outliers at
#'   each step. Default is `0.2`.
#' @param nboot Number of bootstrap replicates to estimate limits more robustly.
#'   Default is `NULL` (no bootstrap, even if **I would recommend to use it**).
#'
#' @return A list with
#' - `x`: the initial vector, whose outliers have been removed,
#' - `lim`: lower and upper limits for outlier removal,
#' - `all_lim`: all bootstrap replicates for `lim` (if `nboot` not `NULL`).
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(1000)
#' str(hist_out(x))
#'
#' # Easy to separate
#' x2 <- c(x, rnorm(50, mean = 7))
#' hist(x2, breaks = nclass.scottRob)
#' str(hist_out(x2))
#'
#' # More difficult to separate
#' x3 <- c(x, rnorm(50, mean = 6))
#' hist(x3, breaks = nclass.scottRob)
#' str(hist_out(x3))
#' str(hist_out(x3, nboot = 999))
#'
hist_out <- function(x, breaks = nclass.scottRob, pmax_out = 0.2, nboot = NULL) {

  if (!is.null(nboot)) {
    boot_repl <- replicate(nboot, {
      x_boot <- sample(x, replace = TRUE)
      hist_out(x_boot, breaks = breaks, pmax_out = pmax_out, nboot = NULL)$lim
    })
    lim <- apply(boot_repl, 1, median)
    return(list(x = x[lim[1] < x & x < lim[2]], lim = lim, all_lim = boot_repl))
  }

  q_inf <- -Inf
  q_sup <- +Inf
  repeat {
    # histogram without outliers
    x_no_out <- x[q_inf < x & x < q_sup]
    h <- graphics::hist(x_no_out, breaks = breaks, plot = FALSE)
    # find empty regions
    ind_wh0 <- which(h$counts == 0)
    val_wh0 <- h$mids[ind_wh0]
    # on the right + define new threshold to refine histogram
    ind_sup_zero <- ind_wh0[val_wh0 > quantile(x_no_out, 1 - pmax_out)]
    if (cont_sup <- (length(ind_sup_zero) > 0))
      q_sup <- h$mids[min(ind_sup_zero)]
    # on the left + define new threshold to refine histogram
    ind_inf_zero <- ind_wh0[val_wh0 < quantile(x_no_out, pmax_out)]
    if (cont_inf <- (length(ind_inf_zero) > 0))
      q_inf <- h$mids[max(ind_inf_zero)]
    if (!cont_sup && !cont_inf) break
  }

  list(x = x_no_out, lim = c(q_inf, q_sup))
}

################################################################################

#' Number of spikes in PCA
#'
#' Estimate the number of distant spikes based on the histogram of eigenvalues.
#'
#' @param eigval Eigenvalues (squared singular values).
#' @param nboot Number of bootstrap replicates to estimate limits more robustly.
#'   Default is `100`.
#' @inheritParams hist_out
#'
#' @return The estimated number of distant spikes.
#'
#' @export
#'
#' @examples
#' N <- 400; M <- 2000; K <- 8
#' U <- matrix(0, N, K); U[] <- rnorm(length(U))
#' V <- matrix(0, M, K); V[] <- rnorm(length(V))
#' # X = U V^T + E
#' X <- tcrossprod(U, V) + 15 * rnorm(N * M)
#' pca <- prcomp(X)
#' eigval <- pca$sdev^2
#' plot(head(eigval, -1), log = "xy", pch = 20)
#' pca_nspike(eigval)
#'
pca_nspike <- function(eigval, breaks = "FD", nboot = 100) {
  lim_up <- hist_out(eigval, breaks = "FD", nboot = nboot)$lim[2]
  sum(eigval > lim_up)
}

################################################################################
