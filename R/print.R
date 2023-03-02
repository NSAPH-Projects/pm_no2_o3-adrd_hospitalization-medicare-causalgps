#' Title
#'
#' @param x TBD
#' @param ... TBD
#'
#' @return
#' value
#' @export
#'
print.path_holder <- function(x, ...){
  x <- unclass(x)
  tmp_names <- names(x)
  for (nm in tmp_names){
    cat(paste("\n", nm, " : ", x[[nm]]))
  }
}
