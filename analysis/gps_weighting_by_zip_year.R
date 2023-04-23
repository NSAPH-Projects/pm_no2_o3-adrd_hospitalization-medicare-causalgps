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

# parameters for this computing job; user should set
n_cores <- 1 # 48 is max of fasse partition, 64 is max of fasse_bigmem partition
n_gb <- 32 # 184 is max of fasse partition, 499 is max of fasse_bigmem partition
find_best_cov_bal_attempt <- F # user should set this variable; true means run for loop over several attempts to find attempt with best covariate balance
save_best_attempt_cov_bal <- F # user should set this variable; true means save covariate balance as csv and plot

if (find_best_cov_bal_attempt){
  n_attempts <- 30 # user should set this; number of attempts this script will try to model the GPS
  n_total_attempts <- 30 # user can set this to a number larger than n_attempts if some attempts have already been tried; to be printed on cov bal plot
  
  if (n_attempts < n_total_attempts){
    modifications <- paste0(n_attempts, "more_attempts") # to be used in names of output files, to record how you're tuning the models
  } else{
    modifications <- paste0(n_attempts, "attempts") # to be used in names of output files, to record how you're tuning the models
  }
} else{
  n_attempts <- 1
  
  if (exposure_name == "pm25"){
    best_maxAC_attempt <- 3
  } else if (exposure_name == "no2"){
    best_maxAC_attempt <- 20
  } else if (exposure_name == "ozone_summer"){
    best_maxAC_attempt <- 27
  } else message("'exposure_name' must be 'pm25', 'no2', or 'ozone_summer'")
  
  modifications <- paste0("attempt", best_maxAC_attempt) # to be used in names of output files, to record how you're tuning the models
}


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
  
  # plot best covariate balance
  weighted_cov_bal_plot <- ggplot(best_maxAC_cov_bal, aes(x = Covariate, y = Absolute_Correlation, color = Dataset, group = Dataset)) +
    geom_point() +
    geom_line() +
    ylab(paste("Absolute Correlation with", exposure_name)) +
    ggtitle(paste0(format(nrow(zip_year_data_with_strata), scientific = F, big.mark = ','), " units of analysis (Attempt #", best_maxAC_attempt, " of ", n_total_attempts, ")")) +
    theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5))
  
  # # save image of best covariate balance plot
  # ggsave(paste0(dir_results, "covariate_balance/",
  #               exposure_name, "/",
  #               "weighting/",
  #               modifications, "/",
  #               nrow(zip_year_data_with_strata), "rows.png"), weighted_cov_bal_plot)
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
  
  # plot best covariate balance
  weighted_cov_bal_plot <- ggplot(best_maxAC_cov_bal, aes(x = Covariate, y = AbsoluteCorrelation, color = Dataset, group = Dataset)) +
    geom_point() +
    geom_line() +
    ylab(paste("Absolute Correlation with", exposure_name)) +
    ggtitle(paste0(format(unique(best_maxAC_cov_bal$SampleSize), scientific = F, big.mark = ','), " units of analysis (Attempt #", best_maxAC_attempt, ")")) +
    theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5))
}

# # run parametric outcome model
# cl <- parallel::makeCluster(n_cores, type = "PSOCK")
# bam_exposure_only <- bam(formula_expos_only,
#                          data = best_weighted_pseudopop,
#                          offset = log(n_persons * n_years),
#                          family = poisson(link = "log"),
#                          weights = capped_stabilized_ipw,
#                          samfrac = 0.05,
#                          chunk.size = 5000,
#                          control = gam.control(trace = TRUE),
#                          nthreads = n_cores,
#                          cluster = cl)
# parallel::stopCluster(cl)

# # get parametric results of interest (coefficient for w)
# coef <- summary(bam_exposure_only)$p.coeff["w"] # alternatively, summary(bam_exposure_only)$p.table["w", "Estimate"]
# coef_se <- summary(bam_exposure_only)$se["w"] # alternatively, summary(bam_exposure_only)$p.table["w", "Std. Error"]
# 
# # save parametric result in specific folder
# parametric_result <- data.table(exposure = exposure_name,
#                                 method = "weighting",
#                                 coefficient = coef,
#                                 se_unadjusted = coef_se)
# fwrite(parametric_result,
#        paste0(dir_results, "parametric_results/",
#               exposure_name, "/",
#               "weighting/",
#               modifications, "/",
#               "parametric_result.csv"))
# 
# # save parametric result in existing table of all parametric results
# parametric_results_table <- fread(paste0(dir_results, "parametric_results/parametric_results_table.csv"))
# parametric_results_table[exposure == exposure_name &
#                            method == "weighting", `:=`(coefficient = coef,
#                                                       se_unadjusted = coef_se)]
# fwrite(parametric_results_table,
#        paste0(dir_results, "parametric_results/parametric_results_table.csv"))

# run semiparametric (thin-plate spline) outcome model
cl <- parallel::makeCluster(n_cores, type = "PSOCK")
bam_exposure_only <- bam(formula_expos_only_smooth,
                         data = best_weighted_pseudopop,
                         offset = log(n_persons * n_years),
                         family = poisson(link = "log"),
                         weights = capped_stabilized_ipw,
                         samfrac = 0.05, # maybe ask if this is good
                         chunk.size = 5000, # can probably make this bigger to run faster but needs more memory
                         control = gam.control(trace = TRUE),
                         nthreads = n_cores,
                         cluster = cl)
parallel::stopCluster(cl)

# # save semiparametric plot
# png(paste0(dir_results, "semiparametric_results/ERFs/",
#            exposure_name, "/",
#            "weighting/",
#            modifications, "/",
#            "bam_smooth_exposure_only.png"))
# plot(bam_exposure_only, main = paste0("GPS ",
#                                       "weighting",
#                                       ", Smoothed Poisson regression,\nexposure only (",
#                                       exposure_name, ")"))
# dev.off()

# Save semiparametric point estimates

# define exposure points at which to predict the outcome
w_values <- seq(min(zip_year_data$w), max(zip_year_data$w), length.out = 20)

# first, use lapply and see if faster than sapply. then, makeCluster(nthread, type = "PSOCK" and use parLapply with cluster
predicted_erf_list <- lapply(w_values,
                             predict_erf_at_a_point,
                             spline_obj = bam_exposure_only,
                             df = best_weighted_pseudopop)

# alternatively, try parLapply
library(parallel)
cl <- parallel::makeCluster(n_cores, type = "PSOCK")
predicted_erf_list <- parLapply(cl = cl,
                                X = w_values,
                                fun = predict_erf_at_a_point,
                                spline_obj = bam_exposure_only,
                                df = best_weighted_pseudopop)
parallel::stopCluster(cl)

# if the above finishes running, save predictions
predicted_erf <- data.table(w = w_values,
                            prediction = unlist(predicted_erf_list))