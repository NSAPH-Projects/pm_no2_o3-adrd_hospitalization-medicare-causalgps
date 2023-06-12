library(data.table)

# get directories and classifications of variables
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/git/code/"
source(paste0(dir_code, "constants.R"))

# read in regression results (coefficient for exposure)
coef <- fread(paste0(dir_results, "parametric_results/coef_for_exposure.txt"))

# read in m-out-of-n bootstrap standard errors (SEs) of coefficient for exposure
boot_SE <- fread(paste0(dir_results, "parametric_results/boot_SE.txt"))
boot_SE[, Method := ifelse(Method == "associational", "Poisson Regression",
                           ifelse(Method == "weighting", "GPS Weighting",
                                  ifelse(Method == "matching", "GPS Matching", "Unknown")))]

# merge point estimates and SEs
all_results <- merge(coef, boot_SE, by = c("Exposure", "Method"))

# get each exposure's IQR
zip_exposure_summary <- fread(paste0(dir_results, "exploratory/zip_exposure_summary.csv"))
iqr <- zip_exposure_summary[, .(Exposure, IQR)] # 4.159 micrograms/m^2 for PM2.5, 12.012 ppb for NO2, 9.801 ppb for summer ozone

# calculate point estimate for hazard ratio (HR) per IQR increase in each exposure
all_results[, HR := exp(IQR * coef)]

# calculate 95% confidence interval for coefficient
all_results[, `:=`(coef_CI_lower = coef - 1.96 * boot_SE,
                   coef_CI_upper = coef + 1.96 * boot_SE)]

# calculate 95% confidence interval for hazard ratio (HR)
all_results[, `:=`(HR_CI_lower = exp(IQR * coef_CI_lower),
                   HR_CI_upper = exp(IQR * coef_CI_upper))]

# save results as csv
fwrite(all_results, paste0(dir_results, "parametric_results/all_HR_results.csv"))
