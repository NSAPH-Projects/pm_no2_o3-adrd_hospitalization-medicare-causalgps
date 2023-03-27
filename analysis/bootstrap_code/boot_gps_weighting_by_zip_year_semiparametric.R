##### Setup #####

# devtools::install_github("fasrc/CausalGPS", ref="develop")
library(data.table)
library(fst)
library(CausalGPS)
library(wCorr)
library(mgcv)
library(tidyr)

# directories for data, code, and results
dir_data <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/data/"
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/code/"
dir_results <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/results/"

# set exposure
exposure_name <- "pm25"

# parameters for this computing job
n_cores <- 1 # 48 is max of fasse partition, 64 is max of fasse_bigmem partition
n_attempts <- 30 # number of attempts for GPS weighting

# get m and n values for this m-out-of-n bootstrap
if (exposure_name == "pm25"){
  n_boot <- 30619
  m_boot <- 2964
} else if (exposure_name == "no2"){
  n_boot <- 30921
  m_boot <- 2990
} else if (exposure_name == "ozone_summer"){
  n_boot <- 30314
  m_boot <- 2937
} else message("'exposure_name' must be 'pm25', 'no2', or 'ozone_summer'")

modifications <- paste0("gps_by_zip_year_", n_attempts, "attempts_boot_",
                        m_boot, "zips") # to be used in names of output files, e.g., cov bal summary

# get helpful functions
source(paste0(dir_code, "analysis/helper_functions.R"))

# get m out of n bootstrap sample
boot_sample_number <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

boot_zips <- read.csv(paste0(dir_data, "analysis/",
                             exposure_name, "/",
                             "boot_zips/",
                             m_boot, "zips_200replicates.csv"))
boot_zips <- boot_zips[, boot_sample_number]

zip_year_data <- read_fst(paste0(dir_data, "analysis/", exposure_name, "/zip_year_data_trimmed_0.05_0.95.fst"), as.data.table = T)
zip_year_data_with_strata <- read_fst(paste0(dir_data, "analysis/", exposure_name, "/zip_year_data_with_strata_trimmed_0.05_0.95.fst"), as.data.table = T)
boot_zip_year_data <- zip_year_data[zip %in% boot_zips]
boot_zip_year_data_with_strata <- zip_year_data_with_strata[zip %in% boot_zips]
rm(zip_year_data_with_strata)

# make sure categorical variables are factors
boot_zip_year_data[, `:=`(zip = as.factor(zip),
                          year = as.factor(year))]
boot_zip_year_data_with_strata[, `:=`(zip = as.factor(zip),
                                      year = as.factor(year),
                                      age_grp = as.factor(age_grp),
                                      sex = as.factor(sex),
                                      race = as.factor(race),
                                      dual = as.factor(dual))]

##### GPS Weighting #####

# set up data.table to check covariate balance for each GPS modeling attempt
cov_bal_weighting <- create_cov_bal_data.table(method = "weighting",
                                               attempt_numbers = 1:n_attempts,
                                               zip_year_data = boot_zip_year_data)

for (i in 1:n_attempts){
  cov_bal_weighting <- get_weighted_pseudopop(attempt_number = i,
                                              zip_year_data = boot_zip_year_data,
                                              zip_year_data_with_strata = boot_zip_year_data_with_strata,
                                              cov_bal_data.table = cov_bal_weighting,
                                              return_cov_bal = T)
}
cov_bal_summary <- summarize_cov_bal(cov_bal_data.table = cov_bal_weighting,
                                     exposure_name = exposure_name,
                                     method = "weighting",
                                     modifications = modifications,
                                     save_csv = F)

# find GPS model with best covariate balance
best_maxAC_attempt <- cov_bal_summary$Attempt[which.min(cov_bal_summary$maxAC)]
best_weighted_pseudopop <- get_weighted_pseudopop(attempt_number = best_maxAC_attempt,
                                                  zip_year_data = boot_zip_year_data,
                                                  zip_year_data_with_strata = boot_zip_year_data_with_strata,
                                                  cov_bal_data.table = cov_bal_weighting,
                                                  return_cov_bal = F)

# fit semiparametric outcome model
bam_exposure_only <- bam(formula_expos_only_smooth,
                         data = best_weighted_pseudopop,
                         offset = log(n_persons * n_years),
                         family = poisson(link = "log"),
                         weights = capped_stabilized_ipw,
                         samfrac = 0.05,
                         chunk.size = 5000,
                         control = gam.control(trace = TRUE),
                         nthreads = n_cores)

# use fitted model to predict log rate of ADRD event at any hypothetical level of exposure
# averaged across all observations in the pseudopopulation
predict_erf_at_a_point <- function(w){
  data <- best_weighted_pseudopop
  data$w <- w
  return(mean(predict(bam_exposure_only, data)))
}

# result to be bootstrapped: predicted values at hypothetical levels of exposure
w_values <- seq(min(zip_year_data$w), max(zip_year_data$w), length.out = 10)
predicted_erf <- sapply(w_values, predict_erf_at_a_point)
predicted_erf <- data.table(w = w_values,
                            prediction = predicted_erf)
fwrite(predicted_erf,
       paste0(dir_results, "bootstrap/",
              exposure_name, "/",
              "weighting/",
              m_boot, "zips/",
              "semiparametric_results/",
              "replicate_", boot_sample_number, ".csv"))