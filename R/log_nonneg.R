
#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
log_nonneg <- function(x){
  if (min(x) >= 0) return(log(x + 0.001)) else return(x)
}


