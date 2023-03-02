#' Title
#'
#' @param x TBD
#'
#' @return
#' value
#' @export
#'
logit_nonneg <- function(x){
  if (min(x) >= 0 & max(x) <= 1){
    return(log((x + 0.001)/(1 - x + 0.001)))
  } else{
    return(x)
  }
}
