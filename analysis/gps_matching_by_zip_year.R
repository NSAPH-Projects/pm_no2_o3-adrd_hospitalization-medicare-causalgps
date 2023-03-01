##### Setup #####

# devtools::install_github("fasrc/CausalGPS", ref="develop")
library(data.table)
library(fst)
library(CausalGPS)
library(wCorr)
library(mgcv)
library(ggplot2)

# directories for data, code, and results
dir_data <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/data/"
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/code/"
dir_results <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/results/"

# set exposure
exposure_name <- "pm25"

# parameters for this computing job
n_cores <- 16 # 48 is max of fasse partition, 64 js max of fasse_bigmem partition
n_gb <- 64 # 184 is max of fasse partition, 499 is max of fasse_bigmem partition
find_best_cov_bal_attempt <- T # user should set this variable
matching_caliper <- 1.5 ## to do: play with this (0.25, 0.5, 1, 1.5). goal is to get large enough ESS

if (find_best_cov_bal_attempt){
  n_attempts <- 10 # user should set this; number of attempts this script will try to model the GPS
  n_total_attempts <- 10 # user can set this to a number larger than n_attempts if some attempts have already been tried; to be printed on cov bal plot
  
  if (n_attempts < n_total_attempts){
    modifications <- paste0(exposure_name, "_only_gps_by_zip_year_trimmed_expos_and_gps_caliper", matching_caliper, "_", n_attempts, "more_attempts") # to be used in names of output files, to record how you're tuning the models
  } else{
    modifications <- paste0(exposure_name, "_only_gps_by_zip_year_trimmed_expos_and_gps_caliper", matching_caliper, "_", n_attempts, "attempts") # to be used in names of output files, to record how you're tuning the models
  }
} else{
  n_attempts <- 1
  best_maxAC_attempt <- 1 # user should set this to the attempt # to be used (for the seed)
  modifications <- paste0(exposure_name, "_only_gps_by_zip_year_trimmed_expos_and_gps_caliper", matching_caliper, "_attempt", best_maxAC_attempt) # to be used in names of output files, to record how you're tuning the models
}

# get data and helpful functions
source(paste0(dir_code, "analysis/helper_functions.R"))
zip_year_data <- read_fst(paste0(dir_data, "analysis/", exposure_name, "_zip_year_data_trimmed_0.05_0.95.fst"),
                          as.data.table = T)
zip_year_data_with_strata <- read_fst(paste0(dir_data, "analysis/", exposure_name, "_zip_year_data_with_strata_trimmed_0.05_0.95.fst"),
                                      as.data.table = T)


##### GPS Matching #####

# get columns from full data that are useful for matching and add a label for stratum
data_for_matching <- copy(zip_year_data_with_strata)
data_for_matching[, stratum := .GRP, by = strata_vars]
setorder(data_for_matching, stratum) # to do: consider if this line is necessary

# set up data.table to check covariate balance for each GPS modeling attempt
if (find_best_cov_bal_attempt){
  n_attempts_already_tried <- n_total_attempts - n_attempts # greater than 0 if user already ran some attempts
  cov_bal_matching <- create_cov_bal_data.table(method = "matching",
                                                 attempt_numbers = (1 + n_attempts_already_tried):n_total_attempts)
} else{
  cov_bal_matching <- create_cov_bal_data.table(method = "matching",
                                                attempt_numbers = best_maxAC_attempt)
}

# function to match within strata; first estimate GPS then match within each stratum
get_matched_pseudopop <- function(attempt_number,
                                  cov_bal_data.table,
                                  return_cov_bal = T,
                                  return_pseudopop = F){
  
  # create log file to see internal processes of CausalGPS
  set_logger(logger_file_path = paste0(dir_code, "analysis/CausalGPS_logs/CausalGPS_", Sys.Date(), "_estimateGpsForMatching_AttemptNumber", attempt_number, "_", modifications, "_", nrow(zip_year_data), "rows_", n_cores, "cores_", n_gb, "gb.log"),
             logger_level = "TRACE")
  
  # set seed according to attempt number
  set.seed(attempt_number)
  
  # estimate GPS
  temp_zip_year_with_gps_obj <- estimate_gps(Y = 0, # fake Y variable since our outcomes are not at the zip-year level; not used in estimate_gps
                                             w = zip_year_data$w,
                                             c = subset(zip_year_data, select = c("year", zip_var_names)),
                                             gps_model = "parametric", # i.e., w=f(x)+epsilon, f(x) estimated by xgboost and epsilon is normal
                                             internal_use = T,
                                             params = list(xgb_nrounds = seq(10, 50),
                                                           xgb_eta = seq(0.1, 0.4, 0.01)),
                                             sl_lib = c("m_xgboost"),
                                             nthread = n_cores)
  
  # create a temporary dataset storing all of CausalGPS's internal parameters, to be expanded from ZIP-years to units of analysis (merging with strata by ZIP, year)
  temp_zip_year_with_gps_dataset_plus_params <- temp_zip_year_with_gps_obj$dataset
  temp_zip_year_with_gps_dataset_plus_params$e_gps_pred <- temp_zip_year_with_gps_obj$e_gps_pred
  temp_zip_year_with_gps_dataset_plus_params$w_resid <- temp_zip_year_with_gps_obj$w_resid
  temp_zip_year_with_gps_dataset_plus_params$zip <- zip_year_data$zip
  
  # trim GPS, keeping observations in middle 95% of GPS values
  gps_outer_quantiles <- quantile(temp_zip_year_with_gps_dataset_plus_params$gps, c(0.025, 0.975))
  temp_zip_year_with_gps_dataset_plus_params <- temp_zip_year_with_gps_dataset_plus_params[gps >= gps_outer_quantiles[1] &
                                                                                             gps <= gps_outer_quantiles[2]]
  
  # apply estimated GPS value to all strata within each remaining/untrimmed ZIP-year (merge will only keep rows in both data.tables)
  # it's important to subset only relevant columns in each data.table:
  # 1) "Y" vector in temp_zip_year_with_gps_dataset_plus_params is fake
  # 2) ZIP-level covariates won't be used in matching or and won't be valid after matching
  # 3) it's okay to include variables for the outcome model because they won't be used in matching anyway and will be useful afterward
  # note that strata_vars includes year
  temp_zip_year_with_gps_dataset_plus_params <- merge(subset(data_for_matching,
                                                             select = c("zip", "Y", "stratum", strata_vars, "n_persons", "n_years")),
                                                      subset(temp_zip_year_with_gps_dataset_plus_params,
                                                             select = c("zip", "year", "w", "gps", "counter_weight", "e_gps_pred", "w_resid")),
                                                      by = c("zip", "year"))
  setorder(temp_zip_year_with_gps_dataset_plus_params, stratum)
  temp_zip_year_with_gps_dataset_plus_params[, row_index := 1:.N, by = stratum]
  temp_zip_year_with_gps_dataset_plus_params$zip <- NULL # ZIP and cohort are the only free variables for matching, should not be used to match
  
  # split data into a list by stratum
  strata_list <- split(temp_zip_year_with_gps_dataset_plus_params,
                       temp_zip_year_with_gps_dataset_plus_params$stratum)
  
  # save which strata are included in the pseudopopulation (since it is possible that some strata may have been entirely trimmed by GPS)
  # names(strata_list) <- unique(temp_zip_year_with_gps_dataset_plus_params$stratum)
  
  # match within strata
  set_logger(logger_file_path = paste0(dir_code, "analysis/CausalGPS_logs/CausalGPS_", Sys.Date(), "_matching_by_stratum", modifications, "_", nrow(zip_year_data_with_strata), "rows_", n_cores, "cores_", n_gb, "gb.log"),
             logger_level = "TRACE")
  temp_matched_pseudopop_list <- lapply(strata_list,
                                        match_within_stratum,
                                        e_gps_std_pred = temp_zip_year_with_gps_obj$e_gps_std_pred,
                                        gps_mx = temp_zip_year_with_gps_obj$gps_mx,
                                        w_mx = temp_zip_year_with_gps_obj$w_mx)
  temp_matched_pseudopop <- rbindlist(temp_matched_pseudopop_list)
  
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

# function to match within a single stratum (dataset_plus_params), after GPS has been estimated
# returns a data.table with variable "counter_weight" denoting number of times each observation is matched
match_within_stratum <- function(dataset_plus_params,
                                 e_gps_std_pred,
                                 gps_mx,
                                 w_mx){
  # make cgps_gps object from input ("dataset_plus_params")
  dataset_as_cgps_gps <- list()
  class(dataset_as_cgps_gps) <- "cgps_gps"
  dataset_as_cgps_gps$dataset <- subset(as.data.frame(dataset_plus_params),
                                        select = c("Y", "w", "gps", "counter_weight", "row_index", # used in matching
                                                   strata_vars, "n_persons", "n_years")) # not used in matching but used in outcome model
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
                                    scale = 1) # notice: max_attempt and transformers are not parameters
  
  return(matched_pop)
}

if (find_best_cov_bal_attempt){
  for (i in 1:n_attempts){
    cov_bal_matching <- get_matched_pseudopop(attempt_number = i + n_attempts_already_tried,
                                              cov_bal_data.table = cov_bal_matching,
                                              return_cov_bal = T,
                                              return_pseudopop = F)
  }
  
  # identify GPS model(s) with best covariate balance
  cov_bal_summary <- summarize_cov_bal(cov_bal_data.table = cov_bal_matching,
                                       method = "matching",
                                       save_csv = T)
  best_maxAC_attempt <- cov_bal_summary$Attempt[which.min(cov_bal_summary$maxAC)]
  best_maxAC_cov_bal <- cov_bal_matching[Attempt == best_maxAC_attempt]
  
  # plot covariate balance
  matched_cov_bal_plot <- ggplot(best_maxAC_cov_bal, aes(x = Covariate, y = Absolute_Correlation, color = Dataset, group = Dataset)) +
    geom_point() +
    geom_line() +
    ylab(paste("Absolute Correlation with", exposure_name)) +
    ggtitle(paste0(format(nrow(zip_year_data_with_strata), scientific = F, big.mark = ','), " units of analysis (Attempt #", best_maxAC_attempt, " of ", n_total_attempts, ")")) +
    theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5))
  
  ggsave(paste0(dir_results, "covariate_balance/matched_pop_", nrow(zip_year_data_with_strata), "rows_", modifications, ".png"), matched_cov_bal_plot)
}

# to do: not lists
# regenerate GPS model and matched pseudopopulation with best covariate balance
best_matched_pseudopop <- get_matched_pseudopop(attempt_number = best_maxAC_attempt,
                                                cov_bal_data.table = cov_bal_matching,
                                                return_cov_bal = F,
                                                return_pseudopop = T)

# print summary statistics for pseudopopulation weights
cat("ESS:", ess(best_matched_pseudopop$counter_weight)) # to do: if ESS is small, investigate which observation(s) are being matched so many times and if increasing? or changing caliper helps
cat("Number of observations matched:", sum(best_matched_pseudopop$counter_weight > 0))
cat("Proportion of observations matched:", sum(best_matched_pseudopop$counter_weight > 0) / nrow(best_matched_pseudopop))
cat("Distribution of number of matches per observations:")
summary(best_matched_pseudopop$counter_weight)
quantile(best_matched_pseudopop$counter_weight, c(0, 0.25, 0.5, 0.75, 0.95, 0.99, 0.999, 0.9999, 1))
boxplot(best_matched_pseudopop$counter_weight)
plot(density(best_matched_pseudopop$w,
             weights = best_matched_pseudopop$counter_weight / sum(best_matched_pseudopop$counter_weight)),
     main = "Density of exposure in matched pseudopopulation",
     xlab = exposure_name)
plot(density(best_matched_pseudopop$w),
     main = "Density of exposure in original population",
     xlab = exposure_name)

# run parametric and semiparametric (thin-plate spline) outcome models
parametric_model_summary <- get_outcome_model_summary(pseudopop = best_matched_pseudopop,
                                                      exposure_name = exposure_name,
                                                      method = "matching",
                                                      n_cores = n_cores,
                                                      parametric_or_semiparametric = "parametric",
                                                      save_results = T)

semiparametric_model_summary <- get_outcome_model_summary(pseudopop = best_matched_pseudopop,
                                                          exposure_name = exposure_name,
                                                          method = "matching",
                                                          n_cores = n_cores,
                                                          parametric_or_semiparametric = "semiparametric",
                                                          save_results = T)