##### Setup #####

library(data.table)
library(fst)

# directories for data, code, and results
dir_data <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/data/"
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/code/"
dir_results <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/results/"

# get data and helpful functions
source(paste0(dir_code, "analysis/helper_functions.R"))
zip_year_data <- read.fst(paste0(dir_data, "analysis/pm_zip_year_data_trimmed_1_99.fst"))

# set up m out of n bootstrap
set.seed(186)
n_boot_iter <- 500
n_boot <- length(unique(zip_year_data$zip))
m_boot <- floor(2 * sqrt(n_boot))
# to do: decide whether with or without replacement
boot_zips <- sapply(1:n_boot_iter, function(i) sample(unique(zip_year_data$zip), m_boot, replace = F))
fwrite(boot_zips, paste0(dir_code, "analysis/batch_code/pm_boot_zips.csv"))
# boot_zips_with_replace <- sapply(1:n_boot_iter, function(i) sample(unique(zip_year_data$zip), m_boot, replace = T))
# fwrite(boot_zips_with_replace, paste0(dir_code, "analysis/pm_boot_zips_with_replace.csv"))