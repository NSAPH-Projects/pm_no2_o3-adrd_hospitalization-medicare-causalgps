
#' Title
#'
#' @param df
#'
#' @return
#' @export
#'
#' @examples
explore_indiv_covs <- function(df){
  # Explore distribution of individual-level covariates
  # To Do: add age distribution
  # Note: offset is not taken into account

  cat("\nProportion male\n")
  print(prop.table(table(df$sexM)))
  print(prop.table(table(df$race_cat)))
  cat("\nProportion Medicaid-eligible\n")
  print(prop.table(table(df$any_dual)))
}
