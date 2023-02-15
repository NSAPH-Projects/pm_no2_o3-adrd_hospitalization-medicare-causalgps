## Classify variables in dataset

offset_var_names <- c("n_persons", "n_years")
zip_expos_names <- c("pm25", "no2", "ozone_summer")
zip_quant_var_names <- c("mean_bmi", "smoke_rate", "hispanic", "prop_blk",
                         "PIR", "poverty", "education", "popdensity", "prop_owner_occ",
                         "summer_tmmx", "summer_rmax")
zip_unordered_cat_var_names <- c("region")
indiv_quant_var_names <- NULL
indiv_unordered_cat_var_names <- c("sex", "race", "dual", "age_grp", "year") # to do: age_group is ordered
strata_vars <- c("year", "sex", "race", "dual", "age_grp") # to do: consider if year should be here
zip_var_names <- c(zip_quant_var_names, zip_unordered_cat_var_names)
indiv_var_names <- c(indiv_unordered_cat_var_names, indiv_quant_var_names) # note: for now, using ADRD_age as a quantitative variable (not binned)

## Formulas for outcome models (parametric and semiparametric thin-plate spline)

formula_expos_only <- as.formula(paste(outcome_name, "~", paste(c("w", strata_vars), collapse = "+", sep = "")))
formula_expos_only_smooth <- as.formula(paste(outcome_name, "~", paste(c("s(w, bs = 'ts')", strata_vars), collapse = "+", sep = "")))


## Functions to assess covariate balance

# to do: see if zip can be included or if need more memory or something
create_cov_bal_data.table <- function(method,
                                      n_attempts,
                                      vars_for_cov_bal = c(other_expos_names, zip_quant_var_names, levels(zip_year_data[["region"]]))){
  
  if (method == "weighting") dataset_names <- c("Weighted", "Unweighted")
  else if (method == "matching") dataset_names <- c("Matched", "Unmatched")
  else stop("'method' must be 'weighting' or 'matching'")
  
  cov_bal <- expand.grid(1:n_attempts,
                         vars_for_cov_bal,
                         dataset_names,
                         100,
                         100) # Correlation and Absolute_Correlation columns will be updated; 100 is placeholder (impossible value of correlation)
  colnames(cov_bal) <- c("Attempt", "Covariate", "Dataset", "Correlation", "Absolute_Correlation")
  cov_bal <- as.data.table(cov_bal)
  return(cov_bal)
}

calculate_correlations <- function(cov_bal_data.table,
                                   method,
                                   attempt,
                                   pseudopop){
  
  if (method == "weighting"){
    dataset_names <- c("Weighted", "Unweighted")
    weight_name <- "capped_stabilized_ipw"
  }
  else if (method == "matching"){
    dataset_names <- c("Matched", "Unmatched")
    weight_name <- "counter_weight"
  }
  else stop("'method' must be 'weighting' or 'matching'")
  
  for (unordered_var in zip_unordered_cat_var_names){
    for (level in levels(zip_year_data[[unordered_var]])){
      cov_bal_data.table[Attempt == attempt & Covariate == level & Dataset == dataset_names[1],
                         Correlation := weightedCorr(pseudopop$w,
                                                     pseudopop[[unordered_var]] == level,
                                                     method = "pearson",
                                                     weights = pseudopop[[weight_name]])]
      
      cov_bal_data.table[Attempt == attempt & Covariate == level & Dataset == dataset_names[2],
                         Correlation := cor(pseudopop$w,
                                            pseudopop[[unordered_var]] == level,
                                            method = "pearson")]
    }
  }
  
  for (quant_var in c(other_expos_names, zip_quant_var_names)){
    cov_bal_data.table[Attempt == attempt & Covariate == quant_var & Dataset == dataset_names[1],
                       Correlation := weightedCorr(pseudopop$w,
                                                   pseudopop[[quant_var]],
                                                   method = "Pearson",
                                                   weights = pseudopop[[weight_name]])]
    
    cov_bal_data.table[Attempt == attempt & Covariate == quant_var & Dataset == dataset_names[2],
                       Correlation := cor(pseudopop$w,
                                          pseudopop[[quant_var]],
                                          method = "pearson")]
  }
  
  return(cov_bal_data.table)
}

summarize_cov_bal <- function(cov_bal_data.table,
                              method,
                              save_csv = T){
  
  if (method == "weighting") dataset_name <- "Weighted"
  else if (method == "matching") dataset_name <- "Matched"
  else stop("'method' must be 'weighting' or 'matching'")
  
  cov_bal_data.table[, Absolute_Correlation := abs(Correlation)]
  cov_bal_summary <- cov_bal_data.table[Dataset == dataset_name, .(maxAC = max(Absolute_Correlation),
                                                                meanAC = mean(Absolute_Correlation),
                                                                maxACVariable = Covariate[which.max(Absolute_Correlation)]), by = Attempt]
  if (save_csv){
    write.csv(cov_bal_summary, paste0(dir_results, "covariate_balance/cov_bal_as_csv/weighted_pop_", n_rows, "rows_", modifications, ".csv"))
  }
  return(cov_bal_summary)
}


## Functions for outcome models

get_outcome_model_summary <- function(pseudopop,
                                      method,
                                      n_cores,
                                      parametric_or_semiparametric = "parametric",
                                      save_results = T){
  
  if (method == "weighting") weight_name = "capped_stabilized_ipw"
  if (method == "matching") weight_name = "counter_weight"
  else stop("'method' must be 'weighting' or 'matching'")
  
  if (parametric_or_semiparametric == "parametric") formula <- formula_expos_only
  else if (parametric_or_semiparametric == "semiparametric") formula <- formula_expos_only_smooth
  else stop("parametric_or_semiparametric must be 'parametric' or 'semiparametric'")
  
  cl <- parallel::makeCluster(n_cores, type = "PSOCK")
  bam_exposure_only <- bam(formula,
                           data = pseudopop,
                           offset = log(n_persons * n_years),
                           family = poisson(link = "log"),
                           weights = pseudopop[[weight_name]],
                           samfrac = 0.05,
                           chunk.size = 5000,
                           control = gam.control(trace = TRUE),
                           nthreads = n_cores,
                           cluster = cl)
  parallel::stopCluster(cl)
  
  if (save_results){
    if (parametric_or_semiparametric == "parametric"){
      saveRDS(summary(bam_exposure_only), file = paste0(dir_results, "parametric_results/bam_exposure_only_", method, "_", n_rows, "rows_", modifications, ".rds"))
    } else{
      png(paste0(dir_results, "semiparametric_results/ERFs/bam_smooth_exposure_only_", method, "_", n_rows, "rows_", modifications, ".png"))
      plot(bam_exposure_only, main = paste0("GPS ", method, ", Smoothed Poisson regression,\nexposure only (", exposure_name, ")"))
      dev.off()
      # saveRDS(bam_exposure_only, file = paste0(dir_results, "semiparametric_results/spline_objects/bam_smooth_exposure_only_", method, "_", n_rows, "rows_", modifications, ".rds"))
    }
  }
  return(summary(bam_exposure_only))
}


## Calculate Kish's effective sample size
ess <- function(weights) return(sum(weights)^2 / (sum(weights^2)))

## Explore distribution of ZIP-level covariates
explore_zip_covs <- function(df){
  cat("\nMin of mean_bmi:", min(df$mean_bmi))
  cat("\nMax of mean_bmi:", max(df$mean_bmi))
  cat("\nMin of smoke_rate:", min(df$smoke_rate))
  cat("\nMax of smoke_rate:", max(df$smoke_rate))
  cat("\nMin of prop_blk:", min(df$prop_blk))
  cat("\nMax of prop_blk:", max(df$prop_blk))
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
  cat("\nMin of prop_owner_occ:", min(df$prop_owner_occ))
  cat("\nMax of prop_owner_occ:", max(df$prop_owner_occ))
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
  # cat("\nMean of prop_blk:", mean(df$prop_blk))
  # cat("\nSD of prop_blk:", sd(df$prop_blk))
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
  # cat("\nMean of prop_owner_occ:", mean(df$prop_owner_occ))
  # cat("\nSD of prop_owner_occ:", sd(df$prop_owner_occ))
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