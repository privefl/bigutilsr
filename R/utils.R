################################################################################

#' Sequence, evenly spaced on a logarithmic scale
#'
#' @inheritParams base::seq
#'
#' @examples
#' seq_log(1, 1000, 4)
#' seq_log(1, 100, 5)
#'
#' @export
seq_log <- function(from, to, length.out) {
  exp(
    seq(from = log(from), to = log(to), length.out = length.out)
  )
}

################################################################################
