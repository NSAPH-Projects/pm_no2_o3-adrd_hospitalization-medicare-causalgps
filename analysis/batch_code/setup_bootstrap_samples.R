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
zip_year_data <- read_fst(paste0(dir_data, "analysis/",
                                 exposure_name, "/",
                                 "zip_year_data_trimmed_0.05_0.95.fst"))

## Set up m out of n bootstrap ZIP codes, with replacement, to be saved as csv

# try multiple possible values of m, as sensitivity analysis
# only requirements of m are that m -> infty and m/n -> 0 as n -> infty
n_boot_iter <- 200
n_boot <- length(unique(zip_year_data$zip)) # 31013?
m_boot_values <- rep(NA, 10)
m_boot_values[1] <- floor(n_boot / log(n_boot))
m_boot_values[2:10] <- sapply(2:10, function(i) floor(i * sqrt(n_boot)))

# for each value of m, bootstrap m ZIP codes and save as csv
for (i in 1:length(m_boot_values)){
  set.seed(186)
  m_boot <- m_boot_values[i]
  boot_zips <- sapply(1:n_boot_iter, function(i) sample(unique(zip_year_data$zip), m_boot, replace = T))
  fwrite(boot_zips, paste0(dir_data,
                           "analysis/",
                           exposure_name, "/",
                           "boot_zips/",
                           m_boot, "zips_",
                           n_boot_iter, "replicates.csv"))
}