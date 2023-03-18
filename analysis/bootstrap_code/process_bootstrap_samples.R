##### Setup #####

library(data.table)

# directories for data, code, and results
dir_data <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/data/"
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/code/"
dir_results <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/results/"

# set exposure
exposure_name <- "pm25"

# set number of bootstrap replicates
n_boot_iter <- 30
n_boot <- 30619 # 30,619 for PM2.5; 30,921 for NO2; 30,314 for ozone
m_boot <- 2964 # user should set this to the desired value of m

# get bootstrap results (1 csv per bootstrap replicate) and merge them
# note this is for weighting for now; to do: use for loop for matching and weighting
boot_results <- lapply(1:n_boot_iter, function(i) fread(paste0(dir_results, "bootstrap/",
                                                               exposure_name, "/",
                                                               "weighting/",
                                                               m_boot, "zips/",
                                                               "replicate_", i, ".csv")))
boot_results <- rbindlist(boot_results)
boot_var <- m_boot / n_boot * var(boot_results$coef_for_exposure,
                                  na.rm = T)
boot_sd <- sqrt(boot_var)