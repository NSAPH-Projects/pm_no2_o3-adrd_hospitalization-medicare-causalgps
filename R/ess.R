
#' Title
#'
#' @param weights
#'
#' @return
#' @export
#'
#' @examples
ess <- function(weights) {
  # Calculate Kish's effective sample size
  return(sum(weights)^2 / (sum(weights^2)))
}
