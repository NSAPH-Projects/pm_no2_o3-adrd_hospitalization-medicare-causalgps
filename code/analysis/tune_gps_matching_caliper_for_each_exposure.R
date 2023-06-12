##### Setup #####

# devtools::install_github("fasrc/CausalGPS", ref="develop")
library(data.table)
library(fst)
library(CausalGPS)
library(wCorr)
library(mgcv)
library(ggplot2)

# get directories, classifications of variables, formulas, and functions
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/git/code/"
source(paste0(dir_code, "constants.R"))
source(paste0(dir_code, "analysis/helper_functions.R"))

# set exposure
exposure_name <- "pm25" # options: "pm25", "no2", or "ozone_summer"

# parameters for this computing job
n_cores <- 4 # 48 is max of fasse partition, 64 is max of fasse_bigmem partition
n_gb <- 48 # 184 is max of fasse partition, 499 is max of fasse_bigmem partition

# get data
zip_year_data <- read_fst(paste0(dir_data, "analysis/",
                                 exposure_name, "/",
                                 "zip_year_data_trimmed_0.05_0.95.fst"),
                          as.data.table = T)
zip_year_data_with_strata <- read_fst(paste0(dir_data, "analysis/",
                                             exposure_name, "/",
                                             "zip_year_data_with_strata_trimmed_0.05_0.95.fst"),
                                      as.data.table = T)

# make sure categorical variables are factors
zip_year_data[, `:=`(zip = as.factor(zip),
                     year = as.factor(year))]
zip_year_data_with_strata[, `:=`(zip = as.factor(zip),
                                 year = as.factor(year),
                                 age_grp = as.factor(age_grp),
                                 sex = as.factor(sex),
                                 race = as.factor(race),
                                 dual = as.factor(dual))]

# set up 5 matching calipers to try for each exposure
if (exposure_name == "pm25"){
  matching_calipers_to_try <- c(0.5, 1, 1.5, 2, 2.5) # recall that IQR for PM2.5 is 4.159
} else if (exposure_name == "no2"){
  matching_calipers_to_try <- c(2.5, 3, 3.5, 4, 4.5) # recall that IQR for NO2 is 12.012
} else if (exposure_name == "ozone_summer"){
  matching_calipers_to_try <- c(2.5, 3, 3.5, 4, 4.5) # recall that IQR for summer ozone is 9.801
} else{
  message("'exposure_name' must be 'pm25', 'no2', or 'ozone_summer'")
}

# each caliper will get 30 attempts to fit a GPS model (each attempt is identified by its seed)
for (matching_caliper in matching_calipers_to_try){
  n_attempts <- 30 # user should set this; number of attempts this script will try to model the GPS
  n_total_attempts <- 30 # user can set this to a number larger than n_attempts if some attempts have already been tried; to be printed on cov bal plot
  
  # a variable called "modifications" will be used in the filepath of the results, to record the covariate balance of each attempt
  if (n_attempts < n_total_attempts){
    modifications <- paste0("match_zips_gps_untrimmed_caliper", matching_caliper, "_", n_attempts, "more_attempts")
  } else{
    modifications <- paste0("match_zips_gps_untrimmed_caliper", matching_caliper, "_", n_attempts, "attempts")
  }
  
  
  ##### GPS Matching #####
  
  # set up data.table to check covariate balance for each GPS modeling attempt
  n_attempts_already_tried <- n_total_attempts - n_attempts # greater than 0 if user already ran some attempts
  cov_bal_matching <- create_cov_bal_data.table(method = "matching",
                                                  attempt_numbers = (1 + n_attempts_already_tried):n_total_attempts,
                                                  zip_year_data = zip_year_data)
  
  # fit GPS model and get matched pseudopopulation using multiple different seeds, to find the best one of 30
  for (i in 1:n_attempts){
    cov_bal_matching <- get_matched_pseudopop(dir_code = dir_code,
                                              attempt_number = i + n_attempts_already_tried,
                                              exposure_name = exposure_name,
                                              modifications = modifications,
                                              n_cores = n_cores,
                                              n_gb = n_gb,
                                              cov_bal_data.table = cov_bal_matching,
                                              zip_year_data = zip_year_data,
                                              zip_year_data_with_strata = zip_year_data_with_strata,
                                              return_cov_bal = T,
                                              return_pseudopop = F)
  }
  
  # identify GPS model(s) with best covariate balance
  cov_bal_summary <- summarize_cov_bal(cov_bal_data.table = cov_bal_matching,
                                       exposure_name = exposure_name,
                                       method = "matching",
                                       modifications = modifications,
                                       save_csv = T)
  best_maxAC_attempt <- cov_bal_summary$Attempt[which.min(cov_bal_summary$maxAC)]
  best_maxAC_cov_bal <- cov_bal_matching[cov_bal_matching$Attempt == best_maxAC_attempt, ]
  
  # save best covariate balance as csv
  fwrite(best_maxAC_cov_bal, paste0(dir_results, "covariate_balance/",
                                    exposure_name, "/",
                                    "matching/",
                                    modifications, "/",
                                    "best_cov_bal.csv"))
}