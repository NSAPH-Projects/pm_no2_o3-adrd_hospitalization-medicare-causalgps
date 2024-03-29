##### Setup #####

library(data.table)
library(fst)

# get directories and classifications of variables
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/git/code/"
source(paste0(dir_code, "constants.R"))

# set exposure
exposure_name <- "pm25" # options: "pm25", "no2", "ozone"

# get data
zip_year_data <- read_fst(paste0(dir_data, "analysis/",
                                 exposure_name, "/",
                                 "zip_year_data_trimmed_0.05_0.95.fst"))

## Set up m out of n bootstrap ZIP codes, with replacement, to be saved as csv

# try multiple possible values of m, as sensitivity analysis
# only requirements of m are that m -> infty and m/n -> 0 as n -> infty
n_boot_iter <- 1000 # or, to start out, 200
n_boot <- length(unique(zip_year_data$zip)) # 30,619 for PM2.5; 30,921 for NO2; 30,314 for ozone
m_boot_values <- rep(NA, 6)
m_boot_values[1:5] <- sapply(1:5, function(i) floor(2 * i * sqrt(n_boot)))
m_boot_values[6] <- floor(n_boot / log(n_boot))

# for each value of m, bootstrap m ZIP codes and save as csv
for (i in 1:length(m_boot_values)){
  set.seed(186)
  m_boot <- m_boot_values[i]
  boot_zips <- sapply(1:n_boot_iter, function(j) sample(unique(zip_year_data$zip), m_boot, replace = T))
  fwrite(boot_zips, paste0(dir_data,
                           "analysis/",
                           exposure_name, "/",
                           "boot_zips/",
                           m_boot, "zips_",
                           n_boot_iter, "replicates.csv"))
}