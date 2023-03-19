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
# note this is for weighting for now; to do: use for loop for matching and weighting
predicted_erf <- lapply(1:n_boot_iter, function(i) fread(paste0(dir_results, "bootstrap/",
                                                               exposure_name, "/",
                                                               "weighting/",
                                                               m_boot, "zips/",
                                                               "semiparametric_results/",
                                                               "replicate_", i, ".csv")))
predicted_erf <- rbindlist(predicted_erf)
predicted_erf_SEs <- predicted_erf[, .(robust_SE = sqrt(m_boot / n_boot * var(prediction))),
                                   by = w] # deleted: na.rm = T

# plot spline with point estimate from full data and SE from m-out-of-n bootstrap

# to do. probably read in point estimate

# ggplot(boot.se,
#        aes(x = x, y = pt.est, 
#            ymin = pt.est - 1.96 * boot.se, ymax = pt.est + 1.96 * boot.se)) +
#   geom_ribbon(fill = "grey", alpha = 0.5) +
#   geom_line() +
#   theme_minimal()