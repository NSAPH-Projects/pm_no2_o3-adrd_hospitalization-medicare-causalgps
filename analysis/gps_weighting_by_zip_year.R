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

# get data and helpful functions
source(paste0(dir_code, "analysis/helper_functions.R"))
zip_year_data <- read_fst(paste0(dir_data, "analysis/", exposure_name, "_zip_year_data_trimmed_1_99.fst"), as.data.table = T)
zip_year_data_with_strata <- read_fst(paste0(dir_data, "analysis/", exposure_name, "_zip_year_data_with_strata_trimmed_1_99.fst"), as.data.table = T)

# parameters for this computing job; user should set
n_cores <- 16 # 48 is max of fasse partition, 64 js max of fasse_bigmem partition
n_gb <- 64 # 184 is max of fasse partition, 499 is max of fasse_bigmem partition
find_best_cov_bal_attempt <- T # user should set this variable

if (find_best_cov_bal_attempt){
  n_attempts <- 30 # user should set this; number of attempts this script will try to model the GPS
  n_total_attempts <- 30 # user can set this to a number larger than n_attempts if some attempts have already been tried; to be printed on cov bal plot
  
  if (n_attempts < n_total_attempts){
    modifications <- paste0(exposure_name, "_only_gps_by_zip_year_", n_attempts, "more_attempts") # to be used in names of output files, to record how you're tuning the models
  } else{
    modifications <- paste0(exposure_name, "_only_gps_by_zip_year_", n_attempts, "attempts") # to be used in names of output files, to record how you're tuning the models
  }
} else{
  n_attempts <- 1
  best_maxAC_attempt <- 1 # user should set this to the attempt # to be used (for the seed)
  modifications <- paste0(exposure_name, "_only_gps_by_zip_year_attempt", best_maxAC_attempt) # to be used in names of output files, to record how you're tuning the models
}


##### GPS Weighting #####

# set up data.table to check covariate balance for each GPS modeling attempt
### to do: see if zip can be included or if need more memory or something
cov_bal_weighting <- create_cov_bal_data.table("weighting", n_attempts)

if (find_best_cov_bal_attempt){
  # create log file to see internal processes of CausalGPS
  set_logger(logger_file_path = paste0(dir_code, "analysis/CausalGPS_logs/CausalGPS_", Sys.Date(), "_estimate_gps_for_weighting_", modifications, "_", nrow(zip_year_data), "rows_", n_cores, "cores_", n_gb, "gb.log"),
             logger_level = "TRACE")
  
  for (i in 1:n_attempts){
    cov_bal_weighting <- get_weighted_pseudopop(attempt_number = i,
                                                zip_year_data = zip_year_data,
                                                zip_year_data_with_strata = zip_year_data_with_strata,
                                                cov_bal_data.table = cov_bal_weighting,
                                                return_cov_bal = T)
  }
  
  cov_bal_summary <- summarize_cov_bal(cov_bal_data.table = cov_bal_weighting,
                                       method = "weighting",
                                       save_csv = T)
  
  # find GPS model with best covariate balance
  best_maxAC_attempt <- cov_bal_summary$Attempt[which.min(cov_bal_summary$maxAC)]
  best_maxAC_cov_bal <- cov_bal_weighting[Attempt == best_maxAC_attempt]
  
  # plot covariate balance
  weighted_cov_bal_plot <- ggplot(best_maxAC_cov_bal, aes(x = Covariate, y = Absolute_Correlation, color = Dataset, group = Dataset)) +
    geom_point() +
    geom_line() +
    ylab(paste("Absolute Correlation with", exposure_name)) +
    ggtitle(paste0(format(nrow(zip_year_data_with_strata), scientific = F, big.mark = ','), " units of analysis (Attempt #", best_maxAC_attempt, " of ", n_total_attempts, ")")) +
    theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5))
  
  ggsave(paste0(dir_results, "covariate_balance/weighted_pop_", nrow(zip_year_data_with_strata), "rows_", modifications, ".png"), weighted_cov_bal_plot)
}

# get pseudopopulation
best_weighted_pseudopop <- get_weighted_pseudopop(attempt_number = best_maxAC_attempt,
                                                  zip_year_data = zip_year_data,
                                                  zip_year_data_with_strata = zip_year_data_with_strata,
                                                  cov_bal_data.table = cov_bal_weighting,
                                                  return_cov_bal = F)

# print summary statistics for pseudopopulation weights
# to do: save in txt file
ess(best_weighted_pseudopop$capped_stabilized_ipw) # for PM2.5 attempt #121, which is best out of 200, ESS is 9,427,355

# run parametric and semiparametric (thin-plate spline) outcome models
parametric_model_summary <- get_outcome_model_summary(pseudopop = best_weighted_pseudopop,
                                                      exposure_name = exposure_name,
                                                      method = "weighting",
                                                      n_cores = n_cores,
                                                      parametric_or_semiparametric = "parametric",
                                                      save_results = T)

semiparametric_model_summary <- get_outcome_model_summary(pseudopop = best_weighted_pseudopop,
                                                          exposure_name = exposure_name,
                                                          method = "weighting",
                                                          n_cores = n_cores,
                                                          parametric_or_semiparametric = "semiparametric",
                                                          save_results = T)