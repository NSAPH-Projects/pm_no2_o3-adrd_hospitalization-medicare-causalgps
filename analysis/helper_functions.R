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

# outcome variable for this analysis
outcome_name <- "n_hosp"


## Formulas for outcome models (parametric and semiparametric thin-plate spline)

formula_expos_only <- as.formula(paste("Y ~", paste(c("w", strata_vars), collapse = "+", sep = "")))
formula_expos_only_smooth <- as.formula(paste("Y ~", paste(c("s(w, bs = 'ts')", strata_vars), collapse = "+", sep = "")))


## Calculate Kish's effective sample size
ess <- function(weights) return(sum(weights)^2 / (sum(weights^2)))


## Functions to assess covariate balance

# to do: see if zip can be included or if need more memory or something
create_cov_bal_data.table <- function(method,
                                      attempt_numbers){
  
  if (method == "weighting") dataset_names <- c("Weighted", "Unweighted")
  else if (method == "matching") dataset_names <- c("Matched", "Unmatched")
  else stop("'method' must be 'weighting' or 'matching'")
  
  vars_for_cov_bal = c(zip_quant_var_names, levels(zip_year_data[["region"]]))
  
  cov_bal <- expand.grid(attempt_numbers,
                         vars_for_cov_bal,
                         dataset_names,
                         -1,
                         -1,
                         -1) # these -1's are placeholders, will be replaced
  colnames(cov_bal) <- c("Attempt", "Covariate", "Dataset", "Correlation", "Absolute_Correlation", "ESS")
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
  
  for (quant_var in zip_quant_var_names){
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
  
  cov_bal_data.table[Attempt == attempt, Absolute_Correlation := abs(Correlation)]
  cov_bal_data.table[Attempt == attempt & Dataset == dataset_names[1], ESS := ess(pseudopop[[weight_name]])]
  cov_bal_data.table[Attempt == attempt & Dataset == dataset_names[2], ESS := nrow(pseudopop)]
  
  return(cov_bal_data.table)
}

summarize_cov_bal <- function(cov_bal_data.table,
                              method,
                              save_csv = T){
  
  if (method == "weighting") dataset_name <- "Weighted"
  else if (method == "matching") dataset_name <- "Matched"
  else stop("'method' must be 'weighting' or 'matching'")
  
  cov_bal_summary <- cov_bal_data.table[Dataset == dataset_name, .(maxAC = max(Absolute_Correlation),
                                                                meanAC = mean(Absolute_Correlation),
                                                                maxACVariable = Covariate[which.max(Absolute_Correlation)],
                                                                ESS = unique(ESS)),
                                        by = Attempt]
  if (save_csv){
    write.csv(cov_bal_summary, paste0(dir_results, "covariate_balance/cov_bal_as_csv/weighted_pop_", modifications, ".csv"))
  }
  return(cov_bal_summary)
}


## Functions to perform GPS weighting or matching

get_weighted_pseudopop <- function(attempt_number,
                                   zip_year_data,
                                   zip_year_data_with_strata,
                                   cov_bal_data.table,
                                   return_cov_bal = T){
  # set seed according to attempt number
  set.seed(attempt_number)
  
  # estimate GPS
  temp_zip_year_with_gps <- estimate_gps(Y = 0, # fake Y variable since our outcomes are not at the zip-year level; not used in estimate_gps
                                         w = zip_year_data$w,
                                         c = subset(zip_year_data, select = c("year", zip_var_names)),
                                         gps_model = "parametric", # i.e., w=f(x)+epsilon, f(x) estimated by xgboost and epsilon is normal
                                         internal_use = T,
                                         params = list(xgb_nrounds = seq(10, 50),
                                                       xgb_eta = seq(0.1, 0.4, 0.01)),
                                         sl_lib = c("m_xgboost"),
                                         nthread = n_cores)
  temp_zip_year_with_gps <- temp_zip_year_with_gps$dataset
  temp_zip_year_with_gps$zip <- zip_year_data$zip
  
  # stabilize GPS using marginal probability of exposure (modeled normally) and cap extreme weights at 10
  marginal_expos_prob <- dnorm(zip_year_data$w,
                               mean = mean(zip_year_data$w),
                               sd = sd(zip_year_data$w))
  temp_zip_year_with_gps$stabilized_ipw <- marginal_expos_prob / temp_zip_year_with_gps$gps ## check estimate_gps
  temp_zip_year_with_gps$capped_stabilized_ipw <- ifelse(temp_zip_year_with_gps$stabilized_ipw > 10, 10, temp_zip_year_with_gps$stabilized_ipw)
  
  # merge with patient data
  # it's important to subset only to relevant columns in temp_zip_year_with_gps because its "Y" vector is fake, so we want to use the real "Y" vector in zip_year_data_with_strata
  temp_weighted_pseudopop <- merge(zip_year_data_with_strata, subset(temp_zip_year_with_gps,
                                                                        select = c("zip", "year", "capped_stabilized_ipw")),
                                   by = c("zip", "year"))
  
  if (return_cov_bal){
    # calculate correlation between exposure and covariates
    cov_bal_weighting <- calculate_correlations(cov_bal_data.table = cov_bal_data.table,
                                                method = "weighting",
                                                attempt = attempt_number,
                                                pseudopop = temp_weighted_pseudopop)
    return(cov_bal_weighting)
  } else{
    return(temp_weighted_pseudopop)
  }
}


## Functions for outcome models

get_outcome_model_summary <- function(pseudopop,
                                      exposure_name,
                                      method,
                                      n_cores,
                                      parametric_or_semiparametric = "parametric",
                                      save_results = T){
  
  if (parametric_or_semiparametric == "parametric") formula <- formula_expos_only
  else if (parametric_or_semiparametric == "semiparametric") formula <- formula_expos_only_smooth
  else stop("parametric_or_semiparametric must be 'parametric' or 'semiparametric'")
  
  cl <- parallel::makeCluster(n_cores, type = "PSOCK")
  if (method == "weighting"){
    bam_exposure_only <- bam(formula,
                             data = pseudopop,
                             offset = log(n_persons * n_years),
                             family = poisson(link = "log"),
                             weights = capped_stabilized_ipw,
                             samfrac = 0.05,
                             chunk.size = 5000,
                             control = gam.control(trace = TRUE),
                             nthreads = n_cores,
                             cluster = cl)
  } else if (method == "matching"){
    bam_exposure_only <- bam(formula,
                             data = pseudopop,
                             offset = log(n_persons * n_years),
                             family = poisson(link = "log"),
                             weights = counter_weight,
                             samfrac = 0.05,
                             chunk.size = 5000,
                             control = gam.control(trace = TRUE),
                             nthreads = n_cores,
                             cluster = cl)
  } else stop("method must be 'weighting' or 'matching'")
  parallel::stopCluster(cl)
  
  if (save_results){
    if (parametric_or_semiparametric == "parametric"){
      saveRDS(summary(bam_exposure_only), file = paste0(dir_results, "parametric_results/bam_exposure_only_", method, "_", nrow(pseudopop), "rows_", modifications, ".rds"))
    } else{
      png(paste0(dir_results, "semiparametric_results/ERFs/bam_smooth_exposure_only_", method, "_", nrow(pseudopop), "rows_", modifications, ".png"))
      plot(bam_exposure_only, main = paste0("GPS ", method, ", Smoothed Poisson regression,\nexposure only (", exposure_name, ")"))
      dev.off()
      # saveRDS(bam_exposure_only, file = paste0(dir_results, "semiparametric_results/spline_objects/bam_smooth_exposure_only_", method, "_", nrow(pseudopop), "rows_", modifications, ".rds"))
    }
  }
  return(summary(bam_exposure_only))
}