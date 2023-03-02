#' Title
#'
#' @param w TBD
#' @param binary_cov TBD
#'
#' @return
#' value
#'
#' @export
#'
abs_pt_biserial_cor <- function(w, binary_cov){
  return(abs(pt_biserial_cor(w, binary_cov)))
}
