
#' Title
#'
#' @param pseudo_pop
#' @param w
#' @param c_unordered_vars
#' @param ci_appr
#' @param all_cov_names
#' @param title
#'
#' @return
#' @export
#'
#' @examples
all_cov_bal <- function(pseudo_pop, w, c_unordered_vars, ci_appr, all_cov_names, title){
  # Check ZIP-level covariate balance in matched data
  # i.e., absolute correlation for quantitative covariates
  # polyserial correlation for ordered categorical variables
  # mean absolute point-biserial correlation for unordered categorical vars
  # params w and c are the same as what was entered into generate_pseudo_pop()
  cor_val_pseudo <- pseudo_pop$adjusted_corr_results$absolute_corr
  cor_val_orig <- pseudo_pop$original_corr_results$absolute_corr

  weights <- pseudo_pop$pseudo_pop$counter_weight

  # correct abs corr values for unordered categorical variables
  for (unordered_var in colnames(c_unordered_vars)){
    cor_val_pseudo[unordered_var] <- weighted_cor_unordered_var(pseudo_pop$pseudo_pop$w, pseudo_pop$pseudo_pop[[unordered_var]], weights)
    cor_val_orig[unordered_var] <- cor_unordered_var(w, c_unordered_vars[[unordered_var]])
  }

  if (ci_appr == "matching"){
    abs_cor = data.frame(Covariate = all_cov_names,
                         Unmatched = cor_val_orig,
                         Matched = cor_val_pseudo) %>%
      gather(c(Unmatched, Matched), key = 'Dataset', value = 'Absolute Correlation')
  } else if (ci_appr == "weighting"){
    abs_cor = data.frame(Covariate = all_cov_names,
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
