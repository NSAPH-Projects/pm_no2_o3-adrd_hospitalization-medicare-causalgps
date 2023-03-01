
#' Title
#'
#' @param w
#' @param unordered_var
#'
#' @return
#' @export
#'
#' @examples
cor_unordered_var <- function(w, unordered_var){
  # Calculate mean of absolute point-biserial correlation between continuous exposure and each binary indicator for an unordered categorical covariate
  # params: w is the vector of continuous exposure, unordered_var is vector of unordered categorical covariate
  levels <- levels(unordered_var) # assumes unordered_var is already a factor, as it should be to be entered into generate_pseudo_pop()
  binary_indicators <- lapply(levels, function(i) 1*(unordered_var == i))
  abs_cor_pb <- lapply(binary_indicators, abs_pt_biserial_cor, w = w)
  return(mean(unlist(abs_cor_pb)))
}
