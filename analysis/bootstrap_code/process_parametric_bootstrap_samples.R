##### Setup #####

library(data.table)

# directories for data, code, and results
dir_data <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/data/"
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/code/"
dir_results <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/results/"

# set exposure and method for this run
exposure_name <- "pm25"
method <- "weighting"

# set number of bootstrap replicates
n_boot_iter <- 200

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

# get bootstrap results (1 csv per bootstrap replicate) and merge them
boot_results <- lapply(1:n_boot_iter, function(i) fread(paste0(dir_results, "bootstrap/",
                                                               exposure_name, "/",
                                                               method, "/",
                                                               m_boot, "zips/",
                                                               "replicate_", i, ".csv")))
boot_results <- rbindlist(boot_results)
boot_var <- m_boot / n_boot * var(boot_results$coef_for_exposure) # deleted: na.rm = T
boot_sd <- sqrt(boot_var)