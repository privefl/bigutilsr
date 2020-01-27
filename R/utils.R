################################################################################

#' Transform a data frame
#'
#' Transform a data frame into a matrix using one hot encoding.
#'
#' @param df A data frame.
#' @param intercept Whether to have a column with all `1`s. Default is `FALSE`.
#'
#' @return A matrix.
#' @export
#'
#' @examples
#' mat <- as_model_matrix(iris)
#' str(mat)
as_model_matrix <- function(df, intercept = FALSE) {

  if (is.data.frame(df)) {
    if (is.null(names(df))) names(df) <- paste0("V", seq_along(df))
    stats::model.matrix(
      stats::as.formula(`if`(intercept, "~ .", "~ . - 1")),
      data = df
    )
  } else {
    stop2("'%s()' is not implemented for class '%s'.",
          as.character(match.call()[[1]]), class(df))
  }
}

################################################################################
