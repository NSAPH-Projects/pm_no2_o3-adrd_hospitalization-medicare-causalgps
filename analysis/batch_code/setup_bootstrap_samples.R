##### Setup #####

library(data.table)
library(fst)

# directories for data, code, and results
dir_data <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/data/"
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/code/"
dir_results <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/results/"

# set exposure
exposure_name <- "pm25"

# get data and helpful functions
source(paste0(dir_code, "analysis/helper_functions.R"))
zip_year_data <- read_fst(paste0(dir_data, "analysis/", exposure_name, "_zip_year_data_trimmed_1_99.fst"))

# set up m out of n bootstrap, with replacement, and save as csv
set.seed(186)
n_boot_iter <- 200
n_boot <- length(unique(zip_year_data$zip)) # 31013
m_boot <- floor(2 * sqrt(n_boot)) # 352
boot_zips <- sapply(1:n_boot_iter, function(i) sample(unique(zip_year_data$zip), m_boot, replace = T))
fwrite(boot_zips, paste0(dir_code, "analysis/batch_code/", exposure_name, "_boot_",
                         m_boot, "zips_",
                         n_boot_iter, "replicates.csv"))

# set up data.frame to store results and save as csv
# note: the following line sets up a filler "zeroth" row, which should be deleted after bootstrap results are generated
boot_results <- data.table(boot_sample_number = 0,
                           coef_for_exposure = 100)
fwrite(boot_results, paste0(dir_results, "batch_sims/", exposure_name, "_boot_results_",
                            m_boot, "zips_",
                            n_boot_iter, "replicates.csv"))