## Calculate Kish's effective sample size
ess <- function(weights) return(sum(weights)^2 / (sum(weights^2)))


## Functions to assess covariate balance

# to do: see if zip can be included or if need more memory or something
create_cov_bal_data.table <- function(method,
                                      attempt_numbers,
                                      zip_year_data){
  
  if (method == "weighting") dataset_names <- c("Weighted", "Unweighted")
  else if (method == "matching") dataset_names <- c("Matched", "Unmatched")
  else stop("'method' must be 'weighting' or 'matching'")
  
  vars_for_cov_bal = c(zip_quant_var_names, levels(zip_year_data[["region"]]))
  
  cov_bal <- expand.grid(attempt_numbers,
                         vars_for_cov_bal,
                         dataset_names,
                         -1,
                         -1,
                         -1,
                         -1,
                         -1) # these -1's are placeholders, will be replaced
  colnames(cov_bal) <- c("Attempt",
                         "Covariate",
                         "Dataset",
                         "Correlation",
                         "AbsoluteCorrelation",
                         "SampleSize",
                         "SampleSizeIncluded", # if matching method, will be smaller than SampleSize
                         "ESS")
  
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
    for (level in levels(pseudopop[[unordered_var]])){
      cov_bal_data.table[Attempt == attempt & Covariate == level & Dataset == dataset_names[1],
                         Correlation := weightedCorr(pseudopop[["w"]],
                                                     pseudopop[[unordered_var]] == level,
                                                     method = "Pearson",
                                                     weights = pseudopop[[weight_name]])]
      
      cov_bal_data.table[Attempt == attempt & Covariate == level & Dataset == dataset_names[2],
                         Correlation := cor(pseudopop[["w"]],
                                            pseudopop[[unordered_var]] == level,
                                            method = "pearson")]
    }
  }
  
  for (quant_var in zip_quant_var_names){
    cov_bal_data.table[Attempt == attempt & Covariate == quant_var & Dataset == dataset_names[1],
                       Correlation := weightedCorr(pseudopop[["w"]],
                                                   pseudopop[[quant_var]],
                                                   method = "Pearson",
                                                   weights = pseudopop[[weight_name]])]
    
    cov_bal_data.table[Attempt == attempt & Covariate == quant_var & Dataset == dataset_names[2],
                       Correlation := cor(pseudopop[["w"]],
                                          pseudopop[[quant_var]],
                                          method = "pearson")]
  }
  
  cov_bal_data.table[Attempt == attempt,
                     AbsoluteCorrelation := abs(Correlation)]
  cov_bal_data.table[Attempt == attempt,
                     SampleSize := nrow(pseudopop)]
  cov_bal_data.table[Attempt == attempt & Dataset == dataset_names[1],
                     SampleSizeIncluded := sum(pseudopop[[weight_name]] > 0)]
  cov_bal_data.table[Attempt == attempt & Dataset == dataset_names[1],
                     ESS := ess(pseudopop[[weight_name]])]
  cov_bal_data.table[Attempt == attempt & Dataset == dataset_names[2],
                     ESS := nrow(pseudopop)]
  
  return(cov_bal_data.table)
}

summarize_cov_bal <- function(cov_bal_data.table,
                              exposure_name,
                              method,
                              modifications,
                              save_csv = T){
  
  if (method == "weighting") dataset_name <- "Weighted"
  else if (method == "matching") dataset_name <- "Matched"
  else stop("'method' must be 'weighting' or 'matching'")
  
  cov_bal_summary <- cov_bal_data.table[Dataset == dataset_name, .(maxAC = max(AbsoluteCorrelation),
                                                                meanAC = mean(AbsoluteCorrelation),
                                                                maxACVariable = Covariate[which.max(AbsoluteCorrelation)],
                                                                SampleSize = unique(SampleSize),
                                                                SampleSizeIncluded = unique(SampleSizeIncluded),
                                                                ESS = unique(ESS)),
                                        by = Attempt]
  if (save_csv){
    write.csv(cov_bal_summary, paste0(dir_results, "covariate_balance/",
                                      exposure_name, "/",
                                      method, "/",
                                      modifications, "/",
                                      "all_cov_bal.csv"))
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

# function to match ZIP codes within years; first estimate GPS on ZIP-year data then match ZIP codes within each year
get_matched_pseudopop <- function(dir_code,
                                  attempt_number,
                                  exposure_name,
                                  modifications,
                                  n_cores,
                                  n_gb,
                                  cov_bal_data.table,
                                  zip_year_data,
                                  zip_year_data_with_strata,
                                  return_cov_bal = T,
                                  return_pseudopop = F){
  
  # # create log file to see internal processes of CausalGPS
  # set_logger(logger_file_path = paste0(dir_code, "analysis/CausalGPS_logs/",
  #                                      exposure_name, "/",
  #                                      "matching/",
  #                                      modifications, "/",
  #                                      Sys.Date(), "_estimateGpsForMatching_AttemptNumber", attempt_number, "_", nrow(zip_year_data), "rows_", n_cores, "cores_", n_gb, "gb.log"),
  #            logger_level = "TRACE")
  
  # set seed according to attempt number
  set.seed(attempt_number)
  
  # estimate GPS
  temp_zip_year_with_gps_obj <- estimate_gps(Y = 0, # fake Y variable since our outcomes are not at the zip-year level; not used in estimate_gps
                                             w = zip_year_data$w,
                                             c = subset(zip_year_data,
                                                        select = c("year", zip_var_names)),
                                             gps_model = "parametric", # i.e., w=f(x)+epsilon, f(x) estimated by xgboost and epsilon is normal
                                             pred_model = "sl",
                                             internal_use = T,
                                             params = list(xgb_nrounds = seq(10, 50),
                                                           xgb_eta = seq(0.1, 0.4, 0.01)),
                                             sl_lib = c("m_xgboost"),
                                             nthread = n_cores)
  
  # create a temporary dataset storing all of CausalGPS's internal parameters, to be split up by year
  temp_zip_year_with_gps_dataset_plus_params <- as.data.table(temp_zip_year_with_gps_obj$dataset)
  temp_zip_year_with_gps_dataset_plus_params$e_gps_pred <- temp_zip_year_with_gps_obj$e_gps_pred
  temp_zip_year_with_gps_dataset_plus_params$w_resid <- temp_zip_year_with_gps_obj$w_resid
  temp_zip_year_with_gps_dataset_plus_params$zip <- zip_year_data$zip
  
  # split dataset into a list of datasets by year
  setkey(temp_zip_year_with_gps_dataset_plus_params, year)
  temp_zip_year_with_gps_dataset_plus_params[, row_index := 1:.N, by = year]
  zip_list_by_year <- split(temp_zip_year_with_gps_dataset_plus_params,
                            temp_zip_year_with_gps_dataset_plus_params$year)
  
  # # create log file to see internal processes of CausalGPS
  # # Note: user must create these folders and subfolders prior to running this line, or else will get error of "No such file or directory"
  # set_logger(logger_file_path = paste0(dir_code, "analysis/CausalGPS_logs/",
  #                                      exposure_name, "/",
  #                                      "matching/",
  #                                      modifications, "/",
  #                                      Sys.Date(), "_matching_zip_by_year", "_", nrow(temp_zip_year_with_gps_dataset_plus_params), "rows_", n_cores, "cores_", n_gb, "gb.log"),
  #            logger_level = "TRACE")
  
  # match within each year
  temp_matched_pseudopop_list <- lapply(zip_list_by_year,
                                        match_zips_within_year,
                                        e_gps_std_pred = temp_zip_year_with_gps_obj$e_gps_std_pred,
                                        gps_mx = temp_zip_year_with_gps_obj$gps_mx,
                                        w_mx = temp_zip_year_with_gps_obj$w_mx)
  
  # combine each year's pseudopopulation to make 1 pseudopopulation for all years
  temp_matched_pseudopop <- rbindlist(temp_matched_pseudopop_list)
  
  # apply estimated GPS value to all strata within each remaining/untrimmed ZIP-year (merge will only keep rows in both data.tables)
  # it's important which columns come from which data.table:
  # 1) "Y" vector in temp_zip_year_with_gps_dataset_plus_params is fake, so use "Y" from data_for_matching
  # 2) include ZIP-level covariates to measure covariate balance, though they may not be valid after matching
  # 3) include variables for the outcome model ("Y", "w", strata_vars, "n_persons", "n_years", "counter_weight")
  # note that strata_vars includes "year"
  temp_matched_pseudopop <- merge(subset(zip_year_data_with_strata,
                                         select = c("zip", "Y", strata_vars, "n_persons", "n_years")), # note that strata_vars includes "year"
                                  subset(temp_matched_pseudopop,
                                         select = c("zip", "year", "w", "counter_weight", zip_var_names)),
                                  by = c("zip", "year"))
  
  if (return_cov_bal){
    cov_bal_data.table <- calculate_correlations(cov_bal_data.table = cov_bal_data.table,
                                                 method = "matching",
                                                 attempt = attempt_number,
                                                 pseudopop = temp_matched_pseudopop)
    return(cov_bal_data.table)
  }
  
  else if (return_pseudopop) return(temp_matched_pseudopop)
  else stop("User must specify whether to return covariate balance or pseudopopulation")
}

# function to match ZIP codes within a single year (dataset_plus_params), after GPS has been estimated
# returns a data.table with variable "counter_weight" denoting number of times each observation is matched
match_zips_within_year <- function(dataset_plus_params,
                                 e_gps_std_pred,
                                 gps_mx,
                                 w_mx){
  # make cgps_gps object from input ("dataset_plus_params")
  dataset_as_cgps_gps <- list()
  class(dataset_as_cgps_gps) <- "cgps_gps"
  dataset_as_cgps_gps$dataset <- subset(as.data.frame(dataset_plus_params),
                                        select = c("Y", "w", "gps", "counter_weight", "row_index", # used in matching
                                                   "zip", "year", # used to merge with full data
                                                   zip_var_names)) # used to measure covariate balance
  dataset_as_cgps_gps$e_gps_pred <- dataset_plus_params$e_gps_pred
  dataset_as_cgps_gps$w_resid <- dataset_plus_params$w_resid
  
  dataset_as_cgps_gps$e_gps_std_pred <- e_gps_std_pred
  dataset_as_cgps_gps$gps_mx <- gps_mx
  dataset_as_cgps_gps$w_mx <- w_mx
  
  matched_pop <- compile_pseudo_pop(data_obj = dataset_as_cgps_gps,
                                    ci_appr = "matching",
                                    gps_model = "parametric",
                                    bin_seq = NULL,
                                    nthread = n_cores,
                                    optimized_compile = T,
                                    matching_fun = "matching_l1",
                                    covar_bl_method = "absolute",
                                    covar_bl_trs = 0.1,
                                    covar_bl_trs_type = "maximal",
                                    delta_n = matching_caliper,
                                    scale = 1) # notice: max_attempt and transformers are not parameters in compile_pseudo_pop(), though they are parameters another commonly used CausalGPS function generate_pseudo_pop()
  
  return(matched_pop)
}


## Functions for outcome models

get_outcome_model_summary <- function(dir_results,
                                      pseudopop,
                                      exposure_name,
                                      method,
                                      modifications,
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
      saveRDS(summary(bam_exposure_only),
              file = paste0(dir_results, "parametric_results/",
                            exposure_name, "/",
                            method, "/",
                            modifications, "/",
                            "bam_exposure_only.rds"))
    } else{
      png(paste0(dir_results, "semiparametric_results/ERFs/",
                 exposure_name, "/",
                 method, "/",
                 modifications, "/",
                 "bam_smooth_exposure_only.png"))
      plot(bam_exposure_only, main = paste0("GPS ",
                                            method,
                                            ", Smoothed Poisson regression,\nexposure only (",
                                            exposure_name, ")"))
      dev.off()
      # saveRDS(bam_exposure_only, file = paste0(dir_results, "semiparametric_results/spline_objects/bam_smooth_exposure_only_", method, "_", nrow(pseudopop), "rows_", modifications, ".rds"))
    }
  }
  return(0)
}

# function to use fitted spline model to predict log rate of ADRD event at any hypothetical level of exposure, averaged across all observations in the population/pseudopopulation
predict_erf_at_a_point <- function(w, spline_obj, df){
  data <- df
  data$w <- w
  return(mean(predict(spline_obj, data)))
}