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
quant_cov_bal <- function(pseudo_pop, ci_appr, var_names, title){
  # Check ZIP-level covariate balance in matched data: abs correlation for quantitative or ordered categorical variables
  cor_val_pseudo <- pseudo_pop$adjusted_corr_results$absolute_corr[var_names] # remove non-ordinal categorical variables; can include zip_ordered_cat_var_names if exists
  cor_val_orig <- pseudo_pop$original_corr_results$absolute_corr[var_names]

  if (ci_appr == "matching"){
    abs_cor = data.frame(Covariate = var_names,
                         Unmatched = cor_val_orig,
                         Matched = cor_val_pseudo) %>%
      gather(c(Unmatched, Matched), key = 'Dataset', value = 'Absolute Correlation')
  } else if (ci_appr == "weighting"){
    abs_cor = data.frame(Covariate = var_names,
                         Unweighted = cor_val_orig,
                         Weighted = cor_val_pseudo) %>%
      gather(c(Unweighted, Weighted), key = 'Dataset', value = 'Absolute Correlation')
  } else stop("ci_appr must be 'matching' or 'weighting'")

  ggplot(abs_cor, aes(x = Covariate, y = `Absolute Correlation`, color = Dataset, group = Dataset)) +
    geom_point() +
    geom_line() +
    ggtitle(title) +
    theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5))
}
