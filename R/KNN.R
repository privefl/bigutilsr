################################################################################

#' @inherit nabor::knn title description params return
#' @inheritDotParams nabor::knn -data -query -k
#' @param ncores Number of cores to use. Default uses [bigparallelr::nb_cores()].
#'
#' @export
#'
#' @examples
#' knn_parallel(matrix(1:4, 2), k = 2, ncores = 2)
#'
knn_parallel <- function(data, query = data, k, ...,
                         ncores = bigparallelr::nb_cores()) {

  all_knn <- bigparallelr::split_parapply(function(ind, data, query, k, ...) {
    nabor::knn(data, query[ind, , drop = FALSE], k = k, ...)
  }, ind = bigparallelr::rows_along(query), ncores = ncores,
  data = data, query = query, k = k, ... = ...)

  list(nn.idx   = do.call("rbind", lapply(all_knn, function(x) x$nn.idx)),
       nn.dists = do.call("rbind", lapply(all_knn, function(x) x$nn.dists)))
}

################################################################################
