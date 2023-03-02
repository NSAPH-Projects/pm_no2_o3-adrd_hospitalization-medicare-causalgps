#' Title
#'
#' @param pseudo_pop TBD
#'
#' @return
#' value
#'
#' @export
#'
summarize_pseudo_weights <- function(pseudo_pop){
  # For GPS-weighted pseudopopulation, print summary statistics for weights
  weights <- pseudo_pop$pseudo_pop$counter_weight
  cat("Number of observations UNTRIMMED by GPS weighting algorithm:", length(pseudo_pop$pseudo_pop$row_index))
  cat("\nNumber of observations with non-zero weight:", sum(weights > 0), "\n")
  cat("\nSum of weights:", sum(weights), "\n")
  print("\nDistribution of weights\n")
  print(quantile(weights, c(0, 0.25, 0.5, 0.75, 0.95, 0.99, 0.999, 1)))
  cat("\nKish ESS:", ess(weights))
}

