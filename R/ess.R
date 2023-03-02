
#' Title
#'
#' @param weights TBD
#'
#' @return
#' value
#' @export
#'
ess <- function(weights) {
  # Calculate Kish's effective sample size
  return(sum(weights)^2 / (sum(weights^2)))
}
