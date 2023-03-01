#' Title
#'
#' @param pseudo_pop
#'
#' @return
#' @export
#'
#' @examples
summarize_pseudo_counter <- function(pseudo_pop){
  # For GPS-matched pseudopopulation, print summary statistics for counter
  counter <- pseudo_pop$pseudo_pop$counter_weight
  cat("Number of observations UNTRIMMED by GPS matching algorithm:", length(pseudo_pop$pseudo_pop$row_index))
  cat("\nNumber of observations matched:", sum(counter > 0), "\n")
  cat("\nNumber of matches:", sum(counter), "\n")
  print("\nDistribution of number of matches per untrimmed observation\n")
  print(quantile(counter, c(0, 0.25, 0.5, 0.75, 0.95, 0.99, 0.999, 1)))
  cat("\nKish ESS:", ess(counter))
}
