# Classify variables in dataset
offset_var_names <- c("n_persons", "n_years")
zip_quant_var_names <- c("mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome",
                         "medianhousevalue", "PIR", "poverty", "education", "popdensity", "pct_owner_occ",
                         "summer_tmmx", "summer_rmax", "no2", "ozone_summer")
zip_unordered_cat_var_names <- c("region", "ADRD_year")
indiv_quant_var_names <- c("ADRD_age")
indiv_unordered_cat_var_names <- c("sexM", "race_cat", "any_dual")
zip_var_names <- c(zip_quant_var_names, zip_unordered_cat_var_names)
indiv_var_names <- c(indiv_unordered_cat_var_names, indiv_quant_var_names) # note: for now, using ADRD_age as a quantitative variable (not binned)

# Functions to transform variables in generate_pseudo_pop()
log_nonneg <- function(x){
  if (min(x) >= 0) return(log(x + 0.001)) else return(x)
}
logit_nonneg <- function(x){
  if (min(x) >= 0 & max(x) <= 1){
    return(log((x + 0.001)/(1 - x + 0.001)))
  } else{
    return(x)
  }
}

# Generate absolute polyserial correlation for each permutation of possible orders of unordered categorical variable and return absolute mean
# If there are more than 24 permutations of the variable's levels (i.e., there are >5 levels), randomly sample
# Note - polyserial correlation is calculated using polycor::polyserial, which "based on the assumption that the joint distribution of the quantitative variable and a latent continuous variable underlying the ordinal variable is bivariate normal"
# param w is the vector of continuous exposures
cor_unordered_var <- function(w, unordered_var, seed = 42){
  library(polycor)
  library(combinat)
  
  abs_polycor_for_1_order <- function(order){
    var_ordered <- factor(unordered_var, levels = order, ordered = T)
    return(abs(polyserial(w, var_ordered)))
  }
  
  levels <- levels(unordered_var)
  
  if (length(levels) > 5){
    set.seed(seed)
    possible_orders <- lapply(1:100, function (i) sample(levels))
  } else possible_orders <- permn(levels)
  
  correlations <- lapply(possible_orders, abs_polycor_for_1_order)
  return(mean(unlist(correlations)))
}

# Check ZIP-level covariate balance in matched data: abs correlation for quantitative or ordered categorical variables, mean polyserial correlation for unordered categorical vars
# params w and c are the same as what was entered into generate_pseudo_pop()
all_cov_bal <- function(pseudo_pop, w, c_unordered_vars, ci_appr, all_cov_names, title){
  cor_val_pseudo <- pseudo_pop$original_corr_results$absolute_corr
  cor_val_orig <- pseudo_pop$adjusted_corr_results$absolute_corr
  
  # correct abs corr values for unordered categorical variables
  for (unordered_var in colnames(c_unordered_vars)){
    cor_val_pseudo[unordered_var] <- cor_unordered_var(pseudo_pop$pseudo_pop$w, pseudo_pop$pseudo_pop[[unordered_var]])
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

# Check ZIP-level covariate balance in matched data: abs correlation for quantitative or ordered categorical variables
quant_cov_bal <- function(pseudo_pop, ci_appr, var_names, title){
  cor_val_pseudo <- pseudo_pop$original_corr_results$absolute_corr[var_names] # remove non-ordinal categorical variables; can include zip_ordered_cat_var_names if exists
  cor_val_orig <- pseudo_pop$adjusted_corr_results$absolute_corr[var_names]
  
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

# Check ZIP-level covariate balance in matched data via boxplot: unordered categorical variables
# to do: contrast test or something
# to do: consider writing function or for loop for all unordered cat vars, or use Xiao's suggestion about polyserial
cat_cov_bal <- function(pseudo_pop, ci_appr, var_names, title){
  if (ci_appr == "matching") weights = "counter"
  else if (ci_appr == "weighting") weights = "ipw"
  else stop("ci_appr must be 'matching' or 'weighting'")
  
  for (var in var_names){
    ggplot(pseudo_pop$pseudo_pop, aes_string(x = var, y = "w", weight = weights)) +
      geom_boxplot() +
      ggtitle(title)
  }
}

# Calculate Kish's effective sample size
ess <- function(weights) return(sum(weights)^2 / (sum(weights^2)))

# Print summary statistics for pseudopopulation counter or weights
summarize_pseudo_weights <- function(pseudo_pop, ci_appr){
  if (ci_appr == "matching") weights = pseudo_pop$pseudo_pop$counter
  else if (ci_appr == "weighting") weights = pseudo_pop$pseudo_pop$ipw
  else stop("ci_appr must be 'matching' or 'weighting'")
  
  # to do: check if these are equivalent
  cat("Number of observations included in pseudo-population:", length(unique(pseudo_pop$pseudo_pop$row_index)))
  cat("Number of observations used in pseudopopulation:", sum(weights > 0))
  
  quantile(weights, c(0, 0.25, 0.5, 0.75, 0.95, 0.99, 0.999))
  cat("Kish ESS:", ess(weights))
}

# Cap counts or weights if desired
cap_weights <- function(pseudo_pop, ci_appr, nthread, quant_var_names, cat_var_names, title){
  if (ci_appr == "matching") weights = "counter"
  else if (ci_appr == "weighting") weights = "ipw"
  else stop("ci_appr must be 'matching' or 'weighting'")
  
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

