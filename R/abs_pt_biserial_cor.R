#' Title
#'
#' @param w
#' @param binary_cov
#'
#' @return
#' @export
#'
#' @examples
abs_pt_biserial_cor <- function(w, binary_cov){
  return(abs(pt_biserial_cor(w, binary_cov)))
}
