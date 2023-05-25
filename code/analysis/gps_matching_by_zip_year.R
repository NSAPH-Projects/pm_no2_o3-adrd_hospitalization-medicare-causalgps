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
exposure_name <- "no2"

# parameters for this computing job
n_cores <- 4 # 48 is max of fasse partition, 64 is max of fasse_bigmem partition
n_gb <- 32 # 184 is max of fasse partition, 499 is max of fasse_bigmem partition
find_best_cov_bal_attempt <- F # user should set this variable; true means run for loop over several attempts to find attempt with best covariate balance
save_best_attempt_cov_bal <- F # user should set this variable; true means save covariate balance as csv and plot

# get matching caliper (previously tuned to the following)
if (exposure_name == "pm25"){
  matching_caliper <- 2
} else if (exposure_name == "no2"){
  matching_caliper <- 3.5
} else if (exposure_name == "ozone_summer"){
  matching_caliper <- 4.5
} else message("'exposure_name' must be 'pm25', 'no2', or 'ozone_summer'")

if (find_best_cov_bal_attempt){
  n_attempts <- 30 # user should set this; number of attempts this script will try to model the GPS
  n_total_attempts <- 30 # user can set this to a number larger than n_attempts if some attempts have already been tried; to be printed on cov bal plot
  
  if (n_attempts < n_total_attempts){
    modifications <- paste0("match_zips_gps_untrimmed_caliper", matching_caliper, "_", n_attempts, "more_attempts") # to be used in names of output files, to record how you're tuning the models
  } else{
    modifications <- paste0("match_zips_gps_untrimmed_caliper", matching_caliper, "_", n_attempts, "attempts") # to be used in names of output files, to record how you're tuning the models
  }
} else{
  n_attempts <- 1
  
  if (exposure_name == "pm25"){
    best_maxAC_attempt <- 23
  } else if (exposure_name == "no2"){
    best_maxAC_attempt <- 8
  } else if (exposure_name == "ozone_summer"){
    best_maxAC_attempt <- 20
  } else message("'exposure_name' must be 'pm25', 'no2', or 'ozone_summer'")
  
  modifications <- paste0("match_zips_gps_untrimmed_caliper", matching_caliper, "_attempt", best_maxAC_attempt) # to be used in names of output files, to record how you're tuning the models
}

# get data and helpful constants and functions
source(paste0(dir_code, "analysis/constants.R"))
source(paste0(dir_code, "analysis/helper_functions.R"))
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


##### GPS Matching #####

# set up data.table to check covariate balance for each GPS modeling attempt
if (find_best_cov_bal_attempt){
  n_attempts_already_tried <- n_total_attempts - n_attempts # greater than 0 if user already ran some attempts
  cov_bal_matching <- create_cov_bal_data.table(method = "matching",
                                                attempt_numbers = (1 + n_attempts_already_tried):n_total_attempts,
                                                zip_year_data = zip_year_data)
} else{
  cov_bal_matching <- create_cov_bal_data.table(method = "matching",
                                                attempt_numbers = best_maxAC_attempt,
                                                zip_year_data = zip_year_data)
}

# if desired, get matched pseudopopulations using multiple different seeds to find the best one
if (find_best_cov_bal_attempt){
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
  
  # plot best covariate balance
  matched_cov_bal_plot <- ggplot(best_maxAC_cov_bal, aes(x = Covariate, y = AbsoluteCorrelation, color = Dataset, group = Dataset)) +
    geom_point() +
    geom_line() +
    ylab(paste("Absolute Correlation with", exposure_name)) +
    ggtitle(paste0(format(unique(best_maxAC_cov_bal$SampleSize), scientific = F, big.mark = ','), " units of analysis (Attempt #", best_maxAC_attempt, " of ", n_total_attempts, ")")) +
    theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5))
  
  # # save image of best covariate balance plot
  # ggsave(paste0(dir_results, "covariate_balance/",
  #               exposure_name, "/",
  #               "matching/",
  #               modifications, "/",
  #               nrow(zip_year_data_with_strata), "rows.png"), matched_cov_bal_plot)
}

# regenerate GPS model and matched pseudopopulation with best covariate balance
best_matched_pseudopop <- get_matched_pseudopop(dir_code = dir_code,
                                                attempt_number = best_maxAC_attempt,
                                                exposure_name = exposure_name,
                                                modifications = modifications,
                                                n_cores = n_cores,
                                                n_gb = n_gb,
                                                cov_bal_data.table = cov_bal_matching,
                                                zip_year_data = zip_year_data,
                                                zip_year_data_with_strata = zip_year_data_with_strata,
                                                return_cov_bal = F,
                                                return_pseudopop = T)

if (save_best_attempt_cov_bal){
  
  # calculate covariate balance of best attempt
  best_maxAC_cov_bal <- calculate_correlations(cov_bal_data.table = cov_bal_matching,
                                               method = "matching",
                                               attempt = best_maxAC_attempt,
                                               pseudopop = best_matched_pseudopop)
  
  # save best covariate balance as csv
  fwrite(best_maxAC_cov_bal, paste0(dir_results, "covariate_balance/",
                                    exposure_name, "/",
                                    "matching/",
                                    modifications, "/",
                                    "best_cov_bal.csv"))
  
  # plot best covariate balance
  matched_cov_bal_plot <- ggplot(best_maxAC_cov_bal, aes(x = Covariate, y = AbsoluteCorrelation, color = Dataset, group = Dataset)) +
    geom_point() +
    geom_line() +
    ylab(paste("Absolute Correlation with", exposure_name)) +
    ggtitle(paste0(format(unique(best_maxAC_cov_bal$SampleSize), scientific = F, big.mark = ','), " units of analysis (Attempt #", best_maxAC_attempt, ")")) +
    theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5))
}


# run parametric outcome model
cl <- parallel::makeCluster(n_cores, type = "PSOCK")
bam_exposure_only <- bam(formula_expos_only,
                         data = best_matched_pseudopop,
                         offset = log(n_persons * n_years),
                         family = poisson(link = "log"),
                         weights = counter_weight,
                         samfrac = 0.05,
                         chunk.size = 5000,
                         control = gam.control(trace = TRUE),
                         nthreads = n_cores,
                         cluster = cl)
parallel::stopCluster(cl)
cat(paste(exposure_name, "GPS Matching", bam_exposure_only$coefficients["w"], sep = ","),
    sep = "\n",
    file = paste0(dir_results, "parametric_results/coef_for_exposure.txt"),
    append = TRUE)

# run semiparametric (thin-plate spline) outcome model
cl <- parallel::makeCluster(n_cores, type = "PSOCK")
bam_exposure_only <- bam(formula_expos_only_smooth_cr,
                         data = best_matched_pseudopop,
                         offset = log(n_persons * n_years),
                         family = poisson(link = "log"),
                         weights = counter_weight,
                         samfrac = 0.05,
                         chunk.size = 5000,
                         control = gam.control(trace = TRUE),
                         nthreads = n_cores,
                         cluster = cl)
parallel::stopCluster(cl)

# estimate counterfactual for every year-zip-strata, calculate ATE
potential_data <- copy(best_matched_pseudopop[, ..strata_vars])
data_prediction <- 
  rbindlist(lapply(seq(min(best_matched_pseudopop$w), 
                       max(best_matched_pseudopop$w), 
                       length.out = 100), function(pot_exp) {
                         
                         # Get potential data if all had same potential exposure
                         potential_data[, w := pot_exp]
                         return(data.table(name = exposure_name,
                                           w = pot_exp,
                                           ate = mean(predict(bam_exposure_only, newdata = potential_data, type = "response"))))
                       }))
plot(I(1e5*ate)~w,data_prediction, type = 'l')
save(data_prediction, file = paste0(dir_results, exposure_name, "_gpsmatching_smooth.rda"))

# save semiparametric point estimates
w_values <- seq(min(zip_year_data$w), max(zip_year_data$w), length.out = 20)
predicted_erf <- sapply(w_values,
                        predict_erf_at_a_point,
                        spline_obj = bam_exposure_only,
                        df = best_matched_pseudopop)
predicted_erf <- data.table(w = w_values,
                            prediction = predicted_erf)
fwrite(predicted_erf,
       paste0(dir_results, "semiparametric_results/",
              exposure_name, "/",
              "matching/",
              modifications, "/",
              "point_estimates.csv"))

# save semiparametric plot
png(paste0(dir_results, "semiparametric_results/ERFs/",
           exposure_name, "/",
           "matching/",
           modifications, "/",
           "bam_smooth_exposure_only.png"))
plot(bam_exposure_only, main = paste0("GPS ",
                                      "matching",
                                      ", Smoothed Poisson regression,\nexposure only (",
                                      exposure_name, ")"))
dev.off()