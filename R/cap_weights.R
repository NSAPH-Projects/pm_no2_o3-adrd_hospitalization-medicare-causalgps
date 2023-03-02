#' Title
#'
#' @param pseudo_pop TBD
#' @param ci_appr TBD
#' @param nthread TBD
#' @param quant_var_names TBD
#' @param cat_var_names TBD
#' @param title TBD
#'
#' @return
#' value
#' @export
#'
cap_weights <- function(pseudo_pop, ci_appr, nthread, quant_var_names, cat_var_names, title){
  # Cap counts or weights if desired
  weights <- "counter_weight"

  cutoff <- quantile(pseudo_pop$pseudo_pop[[weights]], 0.95)
  pseudo_pop$pseudo_pop[[weights]] <- ifelse(pseudo_pop$pseudo_pop[[weights]] > cutoff, cutoff, pseudo_pop$pseudo_pop[[weights]])
  adjusted_corr_obj <- check_covar_balance(pseudo_pop$pseudo_pop,
                                           ci_appr=ci_appr,
                                           nthread=nthread,
                                           covar_bl_method = "absolute",
                                           covar_bl_trs = 0.1,
                                           covar_bl_trs_type = "maximal",
                                           optimized_compile=T)
  pseudo_pop$adjusted_corr_results <-  adjusted_corr_obj$corr_results

  quant_cov_bal(pseudo_pop, ci_appr, quant_var_names, title)
  cat_cov_bal(pseudo_pop, ci_appr, cat_var_names, title)
}
