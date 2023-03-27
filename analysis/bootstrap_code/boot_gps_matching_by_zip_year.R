##### Setup #####

# devtools::install_github("fasrc/CausalGPS", ref="develop")
library(data.table)
library(fst)
library(CausalGPS)
library(wCorr)
library(mgcv)

# directories for data, code, and results
dir_data <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/data/"
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/code/"
dir_results <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/results/"

# set exposure
exposure_name <- "pm25"

# parameters for this computing job
n_cores <- 1 # 48 is max of fasse partition, 64 is max of fasse_bigmem partition
n_gb <- 8 # to do: check if 16 is enough, or if need 32
n_attempts <- 30 # number of attempts for GPS matching

# get m and n values for this m-out-of-n bootstrap
if (exposure_name == "pm25"){
  m_boot <- 2964
} else if (exposure_name == "no2"){
  m_boot <- 2990
} else if (exposure_name == "ozone_summer"){
  m_boot <- 2937
} else message("'exposure_name' must be 'pm25', 'no2', or 'ozone_summer'")

# use same matching caliper as main analysis
if (exposure_name == "pm25"){
  matching_caliper <- 2
} else if (exposure_name == "no2"){
  matching_caliper <- 3.5
} else if (exposure_name == "ozone_summer"){
  matching_caliper <- 4.5
} else message("'exposure_name' must be 'pm25', 'no2', or 'ozone_summer'")

modifications <- paste0("gps_by_zip_year_", n_attempts, "attempts_boot_",
                        m_boot, "zips") # to be used in names of output files, e.g., cov bal summary

# get data and helpful constants and functions
source(paste0(dir_code, "analysis/constants.R"))
source(paste0(dir_code, "analysis/helper_functions.R"))

# get m out of n bootstrap sample
boot_sample_number <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

boot_zips <- read.csv(paste0(dir_data, "analysis/",
                             exposure_name, "/",
                             "boot_zips/",
                             m_boot, "zips_1000replicates.csv"))
boot_zips <- boot_zips[, boot_sample_number]

zip_year_data <- read_fst(paste0(dir_data, "analysis/", exposure_name, "/zip_year_data_trimmed_0.05_0.95.fst"), as.data.table = T)
zip_year_data_with_strata <- read_fst(paste0(dir_data, "analysis/", exposure_name, "/zip_year_data_with_strata_trimmed_0.05_0.95.fst"), as.data.table = T)
boot_zip_year_data <- zip_year_data[zip %in% boot_zips]
boot_zip_year_data_with_strata <- zip_year_data_with_strata[zip %in% boot_zips]
rm(zip_year_data, zip_year_data_with_strata)

# make sure categorical variables are factors
boot_zip_year_data[, `:=`(zip = as.factor(zip),
                          year = as.factor(year))]
boot_zip_year_data_with_strata[, `:=`(zip = as.factor(zip),
                                      year = as.factor(year),
                                      age_grp = as.factor(age_grp),
                                      sex = as.factor(sex),
                                      race = as.factor(race),
                                      dual = as.factor(dual))]

##### GPS Matching #####

# set up data.table to check covariate balance for each GPS modeling attempt
cov_bal_matching <- create_cov_bal_data.table(method = "matching",
                                               attempt_numbers = 1:n_attempts,
                                               zip_year_data = boot_zip_year_data)

for (i in 1:n_attempts){
  cov_bal_matching <- get_matched_pseudopop(dir_code = dir_code,
                                            attempt_number = i,
                                            exposure_name = exposure_name,
                                            modifications = modifications,
                                            n_cores = n_cores,
                                            n_gb = n_gb,
                                            cov_bal_data.table = cov_bal_matching,
                                            zip_year_data = boot_zip_year_data,
                                            zip_year_data_with_strata = boot_zip_year_data_with_strata,
                                            return_cov_bal = T,
                                            return_pseudopop = F)
}
cov_bal_summary <- summarize_cov_bal(cov_bal_data.table = cov_bal_matching,
                                     exposure_name = exposure_name,
                                     method = "matching",
                                     modifications = modifications,
                                     save_csv = F)

# find GPS model with best covariate balance
best_maxAC_attempt <- cov_bal_summary$Attempt[which.min(cov_bal_summary$maxAC)]
best_matched_pseudopop <- get_matched_pseudopop(dir_code = dir_code,
                                                attempt_number = best_maxAC_attempt,
                                                exposure_name = exposure_name,
                                                modifications = modifications,
                                                n_cores = n_cores,
                                                n_gb = n_gb,
                                                cov_bal_data.table = cov_bal_matching,
                                                zip_year_data = boot_zip_year_data,
                                                zip_year_data_with_strata = boot_zip_year_data_with_strata,
                                                return_cov_bal = F,
                                                return_pseudopop = T)

# run parametric outcome model
bam_exposure_only <- bam(formula_expos_only,
                         data = best_matched_pseudopop,
                         offset = log(n_persons * n_years),
                         family = poisson(link = "log"),
                         weights = counter_weight,
                         samfrac = 0.05,
                         chunk.size = 5000,
                         control = gam.control(trace = TRUE),
                         nthreads = n_cores)

# result to be bootstrapped: coefficient for exposure
coef_for_exposure <- summary(bam_exposure_only)$p.table["w", "Estimate"]
result <- data.table(boot_sample_number = boot_sample_number,
                     coef_for_exposure = coef_for_exposure)
fwrite(result,
       paste0(dir_results, "bootstrap/",
              exposure_name, "/",
              "matching/",
              m_boot, "zips/",
              "replicate_", boot_sample_number, ".csv"))
