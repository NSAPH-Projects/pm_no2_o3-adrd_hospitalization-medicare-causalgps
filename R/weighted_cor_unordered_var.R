#' Title
#'
#' @param w
#' @param unordered_var
#' @param weights
#'
#' @return
#' @export
#'
#' @examples
weighted_cor_unordered_var <- function(w, unordered_var, weights){
  levels <- levels(unordered_var) # assumes unordered_var is already a factor, as it should be to be entered into generate_pseudo_pop()
  binary_indicators <- lapply(levels, function(i) 1*(unordered_var == i))
  weighted_cor <- lapply(binary_indicators, weightedCorr, y = w, method = "Pearson", weights = weights)
  abs_weighted_cor <- lapply(weighted_cor, abs)
  return(mean(unlist(abs_weighted_cor)))
}
