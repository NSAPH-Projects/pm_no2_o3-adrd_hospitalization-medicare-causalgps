#' Title
#'
#' @param pseudo_pop TBD
#' @param ci_appr TBD
#' @param var_names TBD
#' @param title TBD
#'
#' @return
#' value
#' @export
#'
cat_cov_bal_boxplot <- function(pseudo_pop, ci_appr, var_names, title){
  # Check ZIP-level covariate balance in matched data via boxplot: unordered categorical variables
  weights <- "counter_weight"

  for (var in var_names){
    ggplot(pseudo_pop$pseudo_pop, aes_string(x = var, y = "w", weight = weights)) +
      geom_boxplot() +
      ggtitle(title)
  }
}
