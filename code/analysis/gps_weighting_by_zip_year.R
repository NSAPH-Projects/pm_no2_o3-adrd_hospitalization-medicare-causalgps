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

# parameters for this computing job; user should set
n_cores <- 4 # 48 is max of fasse partition, 64 is max of fasse_bigmem partition
n_gb <- 48 # 184 is max of fasse partition, 499 is max of fasse_bigmem partition
find_best_cov_bal_attempt <- F # user should set this variable; true means run for loop over several attempts to find attempt with best covariate balance
save_best_attempt_cov_bal <- F # user should set this variable; true means save covariate balance as csv

if (find_best_cov_bal_attempt){
  n_attempts <- 30 # user should set this; number of attempts this script will try to model the GPS
  n_total_attempts <- 30 # user can set this to a number larger than n_attempts if some attempts have already been tried; to be printed on cov bal plot
  
  if (n_attempts < n_total_attempts){
    modifications <- paste0(n_attempts, "more_attempts") # to be used in names of output files, to record how you're tuning the models
  } else{
    modifications <- paste0(n_attempts, "attempts") # to be used in names of output files, to record how you're tuning the models
  }
} else{
  n_attempts <- 1  # the following are the best attempts (out of 30)
  
  if (exposure_name == "pm25"){
    best_maxAC_attempt <- 3
  } else if (exposure_name == "no2"){
    best_maxAC_attempt <- 20
  } else if (exposure_name == "ozone_summer"){
    best_maxAC_attempt <- 27
  } else{
    message("'exposure_name' must be 'pm25', 'no2', or 'ozone_summer'")
  }
  
  modifications <- paste0("attempt", best_maxAC_attempt) # to be used in names of output files, to record how you're tuning the models
}

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


##### GPS Weighting #####

# set up data.table to check covariate balance for each GPS modeling attempt
if (find_best_cov_bal_attempt){
  n_attempts_already_tried <- n_total_attempts - n_attempts # greater than 0 if user already ran some attempts
  cov_bal_weighting <- create_cov_bal_data.table(method = "weighting",
                                                 attempt_numbers = (1 + n_attempts_already_tried):n_total_attempts,
                                                 zip_year_data = zip_year_data)
} else{
  cov_bal_weighting <- create_cov_bal_data.table(method = "weighting",
                                                 attempt_numbers = best_maxAC_attempt,
                                                 zip_year_data = zip_year_data)
}

if (find_best_cov_bal_attempt){
  
  # create log file to see internal processes of CausalGPS
  set_logger(logger_file_path = paste0(dir_code, "analysis/CausalGPS_logs/",
                                       exposure_name, "/",
                                       "weighting/",
                                       modifications, "/",
                                       Sys.Date(), "_estimate_gps_for_weighting_", nrow(zip_year_data), "rows_", n_cores, "cores_", n_gb, "gb.log"),
             logger_level = "TRACE")
  
  for (i in 1:n_attempts){
    cov_bal_weighting <- get_weighted_pseudopop(attempt_number = i + n_attempts_already_tried,
                                                zip_year_data = zip_year_data,
                                                zip_year_data_with_strata = zip_year_data_with_strata,
                                                cov_bal_data.table = cov_bal_weighting,
                                                return_cov_bal = T)
  }
  
  cov_bal_summary <- summarize_cov_bal(cov_bal_data.table = cov_bal_weighting,
                                       exposure_name = exposure_name,
                                       method = "weighting",
                                       modifications = modifications,
                                       save_csv = T)
  
  # find GPS model with best covariate balance
  best_maxAC_attempt <- cov_bal_summary$Attempt[which.min(cov_bal_summary$maxAC)]
  best_maxAC_cov_bal <- cov_bal_weighting[Attempt == best_maxAC_attempt]
  
  # save best covariate balance as csv
  fwrite(best_maxAC_cov_bal, paste0(dir_results, "covariate_balance/",
                                    exposure_name, "/",
                                    "weighting/",
                                    modifications, "/",
                                    "best_cov_bal.csv"))
}

# regenerate GPS model and weighted pseudopopulation with best covariate balance
best_weighted_pseudopop <- get_weighted_pseudopop(attempt_number = best_maxAC_attempt,
                                                  zip_year_data = zip_year_data,
                                                  zip_year_data_with_strata = zip_year_data_with_strata,
                                                  cov_bal_data.table = cov_bal_weighting,
                                                  return_cov_bal = F)

if (save_best_attempt_cov_bal){
  
  # calculate covariate balance of best attempt
  best_maxAC_cov_bal <- calculate_correlations(cov_bal_data.table = cov_bal_weighting,
                                               method = "weighting",
                                               attempt = best_maxAC_attempt,
                                               pseudopop = best_weighted_pseudopop)
  
  # save best covariate balance as csv
  fwrite(best_maxAC_cov_bal, paste0(dir_results, "covariate_balance/",
                                    exposure_name, "/",
                                    "weighting/",
                                    modifications, "/",
                                    "best_cov_bal.csv"))
}

# run parametric outcome model
cl <- parallel::makeCluster(n_cores, type = "PSOCK")
bam_exposure_only <- bam(formula_expos_only,
                         data = best_weighted_pseudopop,
                         offset = log(n_persons * n_years),
                         family = poisson(link = "log"),
                         weights = capped_stabilized_ipw,
                         samfrac = 0.05,
                         chunk.size = 5000,
                         control = gam.control(trace = TRUE),
                         nthreads = n_cores,
                         cluster = cl)
parallel::stopCluster(cl)
cat(paste(exposure_name, "GPS Weighting", bam_exposure_only$coefficients["w"], sep = ","),
    sep = "\n",
    file = paste0(dir_results, "parametric_results/coef_for_exposure.txt"),
    append = TRUE)


### Sensitivity analysis: thin-plate spline outcome model ###

# fit model
cl <- parallel::makeCluster(n_cores, type = "PSOCK")
bam_exposure_only <- bam(formula_expos_only_smooth_cr,
                         data = best_weighted_pseudopop,
                         offset = log(n_persons * n_years),
                         family = poisson(link = "log"),
                         weights = capped_stabilized_ipw,
                         samfrac = 0.05,
                         chunk.size = 5000, # could make this bigger to run faster but needs more memory
                         control = gam.control(trace = TRUE),
                         nthreads = n_cores,
                         cluster = cl)
parallel::stopCluster(cl)

# estimate counterfactual for every year-zip-strata, calculate ATE
potential_data <- copy(best_weighted_pseudopop[, ..strata_vars])
data_prediction <- 
  rbindlist(lapply(seq(min(best_weighted_pseudopop$w), 
                       max(best_weighted_pseudopop$w), 
                       length.out = 100), function(pot_exp) {
    
    # Get potential data if all had same potential exposure
    potential_data[, w := pot_exp]
    return(data.table(name = exposure_name,
                      w = pot_exp,
                      ate = mean(predict(bam_exposure_only, newdata = potential_data, type = "response"))))
  }))
plot(I(1e5*ate)~w,data_prediction, type = 'l')
exposure_density <- density(zip_year_data$w)
save(data_prediction, exposure_density,
     file = paste0(dir_results, "semiparametric_results/", exposure_name, "_gpsweighting_smooth.rda"))