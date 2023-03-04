##### Setup #####

library(data.table)
library(dplyr)

# directories for data, code, and results
dir_data <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/data/"
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/code/"
dir_results <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/results/"

# set exposure
exposure_name <- "pm25"

# set number of bootstrap replicates
n_boot_iter <- 200
n_boot <- 31013 # user should check that this is the correct number of total ZIP codes
m_boot <- 352 # user should set this to the desired value of m; default: floor(2 * sqrt(n_boot))

# get bootstrap results (1 csv per bootstrap replicate) and merge them
# to do: check if the following lines work
boot_results <- list.files(path = paste0(dir_results, "batch_sims/",
                                         exposure_name, "_boot_results/",
                                         m_boot, "zips")) %>% 
  lapply(fread) %>% 
  bind_rows
# boot_results <- boot_results[2:nrow(boot_results), ] # remove initial filler row
boot_var <- m_boot / n_boot * var(boot_results$coef_for_exposure,
                                  na.rm = T)
boot_sd <- sqrt(boot_var)