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

# set up m out of n bootstrap, with replacement, and save as csv
set.seed(186)
n_boot_iter <- 50
n_boot <- length(unique(zip_year_data$zip))
m_boot <- floor(2 * sqrt(n_boot))
boot_zips <- sapply(1:n_boot_iter, function(i) sample(unique(zip_year_data$zip), m_boot, replace = T))
fwrite(boot_zips, paste0(dir_code, "analysis/batch_code/pm_boot_zips.csv"))

# set up data.frame to store results and save as csv
# note: the following line sets up a filler "zeroth" row, which should be deleted after bootstrap results are generated
boot_results <- data.table(boot_sample_number = 0,
                           coef_for_exposure = 100)
fwrite(boot_results, paste0(dir_results, "batch_sims/boot_results.csv"))

##### After running analyses on all bootstrap samples, compute bootstrap estimate and SE #####

# boot_results <- fread(paste0(dir_results, "batch_sims/boot_results.csv"))
# boot_results <- boot_results[2:nrow(boot_results), ] # remove filler row
# boot_mean <- mean(boot_results$coef_for_exposure)
# boot_var <- m_boot / n_boot * var(boot_results$coef_for_exposure)