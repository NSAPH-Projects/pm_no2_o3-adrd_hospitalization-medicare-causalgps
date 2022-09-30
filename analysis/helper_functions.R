# Classify variables in dataset
offset_var_names <- c("n_persons", "n_years")
zip_quant_var_names <- c("mean_bmi", "smoke_rate", "hispanic", "pct_blk",
                         "PIR", "poverty", "education", "popdensity", "pct_owner_occ",
                         "summer_tmmx", "summer_rmax", "no2", "ozone_summer")
zip_unordered_cat_var_names <- c("region", "ADRD_year")
# indiv_quant_var_names <- c("ADRD_age")
indiv_quant_var_names <- NULL
indiv_unordered_cat_var_names <- c("sexM", "race_cat", "any_dual", "ADRD_age_binned")
zip_var_names <- c(zip_quant_var_names, zip_unordered_cat_var_names)
indiv_var_names <- c(indiv_unordered_cat_var_names, indiv_quant_var_names) # note: for now, using ADRD_age as a quantitative variable (not binned)

# Formulas for Poisson regression, including covariate names


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

# Generate point-biserial correlation between continuous exposure (w) and binary covariate (binary_cov)
# Note - equivalent to pearson correlation
pt_biserial_cor <- function(w, binary_cov){
  mean_gp1 <- mean(w[binary_cov == 1])
  mean_gp0 <- mean(w[binary_cov == 0])
  s_nminus1 <- sd(w)
  n1 <- sum(binary_cov == 1)
  n0 <- sum(binary_cov == 0)
  n <- length(binary_cov) # n = n0 + n1
  
  cor_pb <- (mean_gp1 - mean_gp0) / s_nminus1 * sqrt(n1 / n * n0 / (n-1))
  return(cor_pb)
}

abs_pt_biserial_cor <- function(w, binary_cov){
  return(abs(pt_biserial_cor(w, binary_cov)))
}

# Calculate mean of absolute point-biserial correlation between continuous exposure and each binary indicator for an unordered categorical covariate 
# params: w is the vector of continuous exposure, unordered_var is vector of unordered categorical covariate 
cor_unordered_var <- function(w, unordered_var){
  levels <- levels(unordered_var) # assumes unordered_var is already a factor, as it should be to be entered into generate_pseudo_pop()
  binary_indicators <- lapply(levels, function(i) 1*(unordered_var == i))
  abs_cor_pb <- lapply(binary_indicators, abs_pt_biserial_cor, w = w)
  return(mean(unlist(abs_cor_pb)))
}

weighted_cor_unordered_var <- function(w, unordered_var, weights){
  library(wCorr)
  
  levels <- levels(unordered_var) # assumes unordered_var is already a factor, as it should be to be entered into generate_pseudo_pop()
  binary_indicators <- lapply(levels, function(i) 1*(unordered_var == i))
  weighted_cor <- lapply(binary_indicators, weightedCorr, y = w, method = "Pearson", weights = weights)
  abs_weighted_cor <- lapply(weighted_cor, abs)
  return(mean(unlist(abs_weighted_cor)))
}

# Check ZIP-level covariate balance in matched data
# i.e., absolute correlation for quantitative covariates
# polyserial correlation for ordered categorical variables
# mean absolute point-biserial correlation for unordered categorical vars
# params w and c are the same as what was entered into generate_pseudo_pop()
all_cov_bal <- function(pseudo_pop, w, c_unordered_vars, ci_appr, all_cov_names, title){
  cor_val_pseudo <- pseudo_pop$adjusted_corr_results$absolute_corr
  cor_val_orig <- pseudo_pop$original_corr_results$absolute_corr
  
  if (ci_appr == "matching"){
    weights <- pseudo_pop$pseudo_pop$counter
  } else if (ci_appr == "weighting"){
    weights <- pseudo_pop$pseudo_pop$ipw
  } else stop("ci_appr must be 'matching' or 'weighting'")
  
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
cat_cov_bal_boxplot <- function(pseudo_pop, ci_appr, var_names, title){
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

# Explore distribution of ZIP-level covariates
explore_zip_covs <- function(df){
  cat("\nMin of mean_bmi:", min(df$mean_bmi))
  cat("\nMax of mean_bmi:", max(df$mean_bmi))
  cat("\nMin of smoke_rate:", min(df$smoke_rate))
  cat("\nMax of smoke_rate:", max(df$smoke_rate))
  cat("\nMin of pct_blk:", min(df$pct_blk))
  cat("\nMax of pct_blk:", max(df$pct_blk))
  cat("\nMin of hispanic:", min(df$hispanic))
  cat("\nMax of hispanic:", max(df$hispanic))
  cat("\nMin of education:", min(df$education))
  cat("\nMax of education:", max(df$education))
  cat("\nMin of popdensity:", min(df$popdensity))
  cat("\nMax of popdensity:", max(df$popdensity))
  cat("\nMin of poverty:", min(df$poverty))
  cat("\nMax of poverty:", max(df$poverty))
  cat("\nMin of PIR:", min(df$PIR))
  cat("\nMax of PIR:", max(df$PIR))
  cat("\nMin of pct_owner_occ:", min(df$pct_owner_occ))
  cat("\nMax of pct_owner_occ:", max(df$pct_owner_occ))
  cat("\nMin of summer_tmmx:", min(df$summer_tmmx))
  cat("\nMax of summer_tmmx:", max(df$summer_tmmx))
  cat("\nMin of summer_rmax:", min(df$summer_rmax))
  cat("\nMax of summer_rmax:", max(df$summer_rmax))
  cat("\nMin of no2:", min(df$no2))
  cat("\nMax of no2:", max(df$no2))
  cat("\nMin of ozone_summer:", min(df$ozone_summer))
  cat("\nMax of ozone_summer:", max(df$ozone_summer))
  prop.table(table(df$region))
  
  # cat("\nMean of mean_bmi:", mean(df$mean_bmi))
  # cat("\nSD of mean_bmi:", sd(df$mean_bmi))
  # cat("\nMean of smoke_rate:", mean(df$smoke_rate))
  # cat("\nSD of smoke_rate:", sd(df$smoke_rate))
  # cat("\nMean of pct_blk:", mean(df$pct_blk))
  # cat("\nSD of pct_blk:", sd(df$pct_blk))
  # cat("\nMean of hispanic:", mean(df$hispanic))
  # cat("\nSD of hispanic:", sd(df$hispanic))
  # cat("\nMean of education:", mean(df$education))
  # cat("\nSD of education:", sd(df$education))
  # cat("\nMean of popdensity:", mean(df$popdensity))
  # cat("\nSD of popdensity:", sd(df$popdensity))
  # cat("\nMean of poverty:", mean(df$poverty))
  # cat("\nSD of poverty:", sd(df$poverty))
  # cat("\nMean of medhouseholdincome:", mean(df$medhouseholdincome))
  # cat("\nSD of medhouseholdincome:", sd(df$medhouseholdincome))
  # cat("\nMean of PIR:", mean(df$PIR))
  # cat("\nSD of PIR:", sd(df$PIR))
  # cat("\nMean of pct_owner_occ:", mean(df$pct_owner_occ))
  # cat("\nSD of pct_owner_occ:", sd(df$pct_owner_occ))
}

# Explore distribution of individual-level covariates
# To Do: add age distribution
# Note: offset is not taken into account
explore_indiv_covs <- function(df){
  cat("\nProportion male\n")
  print(prop.table(table(df$sexM)))
  print(prop.table(table(df$race_cat)))
  cat("\nProportion Medicaid-eligible\n")
  print(prop.table(table(df$any_dual)))
}

# For GPS-matched pseudopopulation, print summary statistics for counter
summarize_pseudo_counter <- function(pseudo_pop){
  counter <- pseudo_pop$pseudo_pop$counter
  cat("Number of observations UNTRIMMED by GPS matching algorithm:", length(pseudo_pop$pseudo_pop$row_index))
  cat("\nNumber of observations matched:", sum(counter > 0), "\n")
  cat("\nNumber of matches:", sum(counter), "\n")
  print("\nDistribution of number of matches per untrimmed observation\n")
  print(quantile(counter, c(0, 0.25, 0.5, 0.75, 0.95, 0.99, 0.999, 1)))
  cat("\nKish ESS:", ess(counter))
}

# For GPS-weighted pseudopopulation, print summary statistics for weights
summarize_pseudo_weights <- function(pseudo_pop){
  weights <- pseudo_pop$pseudo_pop$ipw
  cat("Number of observations UNTRIMMED by GPS weighting algorithm:", length(pseudo_pop$pseudo_pop$row_index))
  cat("\nNumber of observations with non-zero weight:", sum(weights > 0), "\n")
  cat("\nSum of weights:", sum(weights), "\n")
  print("\nDistribution of weights\n")
  print(quantile(weights, c(0, 0.25, 0.5, 0.75, 0.95, 0.99, 0.999, 1)))
  cat("\nKish ESS:", ess(weights))
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

