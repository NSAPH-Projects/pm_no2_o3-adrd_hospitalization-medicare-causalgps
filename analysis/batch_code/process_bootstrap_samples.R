##### Setup #####

library(data.table)

# directories for data, code, and results
dir_data <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/data/"
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/code/"
dir_results <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/results/"

# set exposure
exposure_name <- "pm25"

# get original data
zip_year_data <- read_fst(paste0(dir_data, "analysis/", exposure_name, "_zip_year_data_trimmed_1_99.fst"))

# set number of bootstrap replicates
n_boot_iter <- 200
n_boot <- length(unique(zip_year_data$zip)) # 31013
m_boot <- floor(2 * sqrt(n_boot)) # 352

# get bootstrap results
boot_results <- read.csv(paste0(dir_results, "batch_sims/", exposure_name, "_boot_results_",
                                m_boot, "zips_",
                                n_boot_iter, "replicates.csv"),
                         fill = T)
boot_results <- boot_results[2:nrow(boot_results), ] # remove initial filler row
boot_var <- m_boot / n_boot * var(boot_results$coef_for_exposure,
                                  na.rm = T)
boot_sd <- sqrt(boot_var)