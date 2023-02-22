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

# get helpful functions
source(paste0(dir_code, "analysis/helper_functions.R"))

# parameters for this computing job
n_cores <- 1 # 48 is max of fasse partition, 64 is max of fasse_bigmem partition
n_attempts <- 30
n_total_attempts <- n_attempts # user can set this to a number larger than n_attempts if some attempts with different seeds have already been tried
modifications <- paste0("gps_by_zip_year_", n_attempts, "attempts_boot_test") # to be used in names of output files, e.g., cov bal summary

# get m out of n bootstrap sample
boot_sample_number <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

boot_zips <- read.csv(paste0(dir_code, "analysis/batch_code/pm_boot_zips.csv"))
boot_zips <- boot_zips[, boot_sample_number]

zip_year_data <- read.fst(paste0(dir_data, "analysis/pm_zip_year_data_trimmed_1_99.fst"), as.data.table = T)
zip_year_data_with_strata <- read.fst(paste0(dir_data, "analysis/ADRD_agg_lagged_trimmed_1_99.fst"), as.data.table = T)
boot_zip_year_data <- zip_year_data[zip %in% boot_zips]
boot_zip_year_data_with_strata <- zip_year_data_with_strata[zip %in% boot_zips]


##### GPS Weighting #####

# set up data.table to check covariate balance for each GPS modeling attempt
### to do: see if zip can be included or if need more memory or something
cov_bal_weighting <- create_cov_bal_data.table("weighting", n_attempts)

for (i in 1:n_attempts){
  cov_bal_weighting <- get_weighted_pseudopop(i,
                                              boot_zip_year_data,
                                              boot_zip_year_data_with_strata,
                                              cov_bal_weighting,
                                              return_cov_bal = T)
}
cov_bal_summary <- summarize_cov_bal(cov_bal_data.table = cov_bal_weighting,
                                     method = "weighting",
                                     save_csv = F)

# find GPS model with best covariate balance
best_maxAC_attempt <- cov_bal_summary$Attempt[which.min(cov_bal_summary$maxAC)]
best_weighted_pseudopop <- get_weighted_pseudopop(best_maxAC_attempt,
                                                  boot_zip_year_data,
                                                  boot_zip_year_data_with_strata,
                                                  cov_bal_weighting,
                                                  return_cov_bal = F)

# run parametric outcome model
weights <- best_weighted_pseudopop$capped_stabilized_ipw # note: to use the following function, need to have "weights" in global environment; to do: improve this
parametric_model_summary <- get_outcome_model_summary(pseudopop = best_weighted_pseudopop,
                                                      method = "weighting",
                                                      n_cores = n_cores,
                                                      parametric_or_semiparametric = "parametric",
                                                      save_results = F)

# result to be bootstrapped: coefficient for exposure
coef_for_exposure <- parametric_model_summary$p.table["w", "Estimate"]
result <- data.table(boot_sample_number = boot_sample_number,
                     coef_for_exposure = coef_for_exposure)
fwrite(result,
       paste0(dir_results, "batch_sims/boot_results.csv"),
       append = T)


##### After running analyses on all bootstrap samples, compute bootstrap estimate and SE #####

# boot_mean <- mean(boot_results)
# boot_var <- m_boot / n_boot * var(boot_results)