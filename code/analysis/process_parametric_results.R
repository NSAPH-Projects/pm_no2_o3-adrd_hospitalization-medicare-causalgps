library(data.table)

# directories for data, code, and results
dir_data <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/data/"
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/code/"
dir_results <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/results/"

# read in regression results (coefficient for exposure)
coef <- fread(paste0(dir_results, "parametric_results/coef_for_exposure.txt"))

# read in m-out-of-n bootstrap standard errors (SEs) of coefficient for exposure
boot_SE <- fread(paste0(dir_results, "parametric_results/boot_SE.txt"))

# merge point estimates and SEs
all_results <- merge(boef, boot_SE, by = c("Exposure", "Method"))

# calculate point estimate for hazard ratio (HR)
all_results$HR <- 1
all_results[Exposure == "pm25", HR := exp(4.159 * coef)]
all_results[Exposure == "no2", HR := exp(12.012 * coef)]
all_results[Exposure == "ozone_summer", HR := exp(9.801 * coef)]

# calculate 95% confidence interval for coefficient
all_results[, `:=`(coef_CI_lower = coef - 1.96 * boot_SE,
                   coef_CI_upper = coef + 1.96 * boot_SE)]

# calculate 95% confidence interval for hazard ratio (HR)
all_results$HR_CI_lower <- 1
all_results$HR_CI_upper <- 1
all_results[Exposure == "pm25", `:=`(HR_CI_lower = exp(4.159 * coef_CI_lower),
                                     HR_CI_upper = exp(4.159 * coef_CI_upper))]
all_results[Exposure == "no2", `:=`(HR_CI_lower = exp(12.012 * coef_CI_lower),
                                    HR_CI_upper = exp(12.012 * coef_CI_upper))]
all_results[Exposure == "ozone_summer", `:=`(HR_CI_lower = exp(9.801 * coef_CI_lower),
                                             HR_CI_upper = exp(9.801 * coef_CI_upper))]