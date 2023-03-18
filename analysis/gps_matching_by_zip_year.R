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
  n_attempts <- 30 # user should set this; number of attempts this script will try to model the GPS
  n_total_attempts <- 30 # user can set this to a number larger than n_attempts if some attempts have already been tried; to be printed on cov bal plot
  
  if (n_attempts < n_total_attempts){
    modifications <- paste0("match_zips_caliper", matching_caliper, "_", n_attempts, "more_attempts") # to be used in names of output files, to record how you're tuning the models
  } else{
    modifications <- paste0("match_zips_caliper", matching_caliper, "_", n_attempts, "attempts") # to be used in names of output files, to record how you're tuning the models
  }
} else{
  n_attempts <- 1
  best_maxAC_attempt <- 1 # user should set this to the attempt # to be used (for the seed)
  modifications <- paste0("match_zips_caliper", matching_caliper, "_attempt", best_maxAC_attempt) # to be used in names of output files, to record how you're tuning the models
}

# get data and helpful functions
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
  
  # plot covariate balance
  matched_cov_bal_plot <- ggplot(best_maxAC_cov_bal, aes(x = Covariate, y = Absolute_Correlation, color = Dataset, group = Dataset)) +
    geom_point() +
    geom_line() +
    ylab(paste("Absolute Correlation with", exposure_name)) +
    ggtitle(paste0(format(unique(best_maxAC_cov_bal$SampleSize), scientific = F, big.mark = ','), " units of analysis (Attempt #", best_maxAC_attempt, " of ", n_total_attempts, ")")) +
    theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5))
  
  ggsave(paste0(dir_results, "covariate_balance/",
                exposure_name, "/",
                "matching/",
                modifications, "/",
                nrow(zip_year_data_with_strata), "rows.png"), matched_cov_bal_plot)
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

# print summary statistics for pseudopopulation weights
cat("ESS:", ess(best_matched_pseudopop$counter_weight)) # to do: if ESS is small, investigate which observation(s) are being matched so many times and if increasing? or changing caliper helps
cat("Number of observations matched:", sum(best_matched_pseudopop$counter_weight > 0))
cat("Proportion of observations matched:", sum(best_matched_pseudopop$counter_weight > 0) / nrow(best_matched_pseudopop))
cat("Distribution of number of matches per observations:")
summary(best_matched_pseudopop$counter_weight)
quantile(best_matched_pseudopop$counter_weight, c(0, 0.25, 0.5, 0.75, 0.95, 0.99, 0.999, 0.9999, 1))
# boxplot(best_matched_pseudopop$counter_weight)

# print number of observations not trimmed by GPS
cat("Number of observations not trimmed for having extreme GPS:", nrow(best_matched_pseudopop))

# print summary of pseudopopulation exposure
best_matched_pseudopop[counter_weight > 0, .(max_exposure = max(w))]
best_matched_pseudopop[counter_weight > 0, .(min_exposure = min(w))]
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
                                                      modifications = modifications,
                                                      n_cores = n_cores,
                                                      parametric_or_semiparametric = "parametric",
                                                      save_results = T)

semiparametric_model_summary <- get_outcome_model_summary(pseudopop = best_matched_pseudopop,
                                                          exposure_name = exposure_name,
                                                          method = "matching",
                                                          modifications = modifications,
                                                          n_cores = n_cores,
                                                          parametric_or_semiparametric = "semiparametric",
                                                          save_results = T)