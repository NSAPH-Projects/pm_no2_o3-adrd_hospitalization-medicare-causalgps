##### Setup #####

library(data.table)

# get directories and classifications of variables
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/git/code/"
source(paste0(dir_code, "constants.R"))


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

# save results in results folder
cat(paste(exposure_name, method, boot_sd, sep = ","),
    sep = "\n",
    file = paste0(dir_results, "parametric_results/boot_SE.txt"),
    append = TRUE)