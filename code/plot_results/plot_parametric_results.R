library(data.table)
library(ggplot2)

# directories for data, code, and results
dir_data <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/data/"
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/code/"
dir_results <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/results/"

# get helpful constants
source(paste0(dir_code, "analysis/constants.R"))

# # get exposure IQRs (for hazard ratio)
# zip_exposure_summary <- fread(paste0(dir_results, "exploratory/zip_exposure_summary.csv"))

# # get results
# parametric_results_table <- fread(paste0(dir_results, "parametric_results/parametric_results_table.csv"))
# 
# # get 95% confidence intervals for regression coefficient
# parametric_results_table[, `:=`(coef_95CI_lower = coefficient - 1.96 * se_bootstrap,
#                                 coef_95CI_upper = coefficient + 1.96 * se_bootstrap)]
# 
# # convert coefficient to hazard ratio
# parametric_results_table$HR <- 0.1 # placeholder value; a double
# parametric_results_table[exposure == "pm25", HR := exp(zip_exposure_summary[Exposure == "pm25", IQR] * coefficient)]
# parametric_results_table[exposure == "no2", HR := exp(zip_exposure_summary[Exposure == "no2", IQR] * coefficient)]
# parametric_results_table[exposure == "ozone_summer", HR := exp(zip_exposure_summary[Exposure == "ozone_summer", IQR] * coefficient)]
# 
# # get 95% confidence intervals for HRs
# parametric_results_table[exposure == "pm25", `:=`(HR_95CI_lower = exp(zip_exposure_summary[Exposure == "pm25", IQR] * coef_95CI_lower),
#                                                   HR_95CI_upper = exp(zip_exposure_summary[Exposure == "pm25", IQR] * coef_95CI_upper))]
# parametric_results_table[exposure == "no2", `:=`(HR_95CI_lower = exp(zip_exposure_summary[Exposure == "no2", IQR] * coef_95CI_lower),
#                                                  HR_95CI_upper = exp(zip_exposure_summary[Exposure == "no2", IQR] * coef_95CI_upper))]
# parametric_results_table[exposure == "ozone_summer", `:=`(HR_95CI_lower = exp(zip_exposure_summary[Exposure == "ozone_summer", IQR] * coef_95CI_lower),
#                                                           HR_95CI_upper = exp(zip_exposure_summary[Exposure == "ozone_summer", IQR] * coef_95CI_upper))]

# for now, hard code results from google sheet
parametric_results_table <- expand.grid(Exposure = c("PM2.5", "NO2", "Summer Ozone"), # note: important to NOT call this column "exposure_name", or else won't be able to match with results later
                                        Method = c("Poisson Regression", "GPS Weighting", "GPS Matching"),
                                        HR = -1000.1,
                                        HR_95CI_lower = -1000.1,
                                        HR_95CI_upper = -1000.1) # these are placeholder numbers (important that they are doubles)
parametric_results_table <- setorder(parametric_results_table, Exposure)
parametric_results_table$HR <- c(1.152275543,
                                 1.108096989,
                                 1.110449857,
                                 1.144553815,
                                 1.05769926,
                                 1.069683249,
                                 1.056699348,
                                 1.044876482,
                                 1.021095561)
parametric_results_table$HR_95CI_lower <- c(1.144153628,
                                            1.09705365,
                                            1.081006371,
                                            1.13738924,
                                            1.048598491,
                                            1.052149766,
                                            1.049753462,
                                            1.035500191,
                                            1.005404529)
parametric_results_table$HR_95CI_upper <- c(1.160455111,
                                            1.119251495,
                                            1.140695299,
                                            1.151763521,
                                            1.066879016,
                                            1.087508917,
                                            1.063691192,
                                            1.054337674,
                                            1.037031478)

# plot results
ggplot(parametric_results_table, aes(x = Exposure,
                                          y = HR,
                                          color = Method,
                                          shape = Method)) +
  geom_hline(yintercept = 1, linewidth = 1, linetype = 2) +
  geom_point(position=position_dodge(0.5), size = 4) +
  geom_errorbar(aes(ymin = HR_95CI_lower,
                    ymax = HR_95CI_upper),
                position=position_dodge(0.5), width = 0.4,
                linewidth = 1) +
  ylab("Hazard Ratio per IQR Increase in Exposure") +
  theme_minimal(base_size = 24) +
  theme(legend.key.width = unit(1, "cm"), 
        legend.key.height = unit(1, "cm"),
        legend.box.spacing = unit(-1, "cm"),
        legend.position = "bottom",
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(color = "", shape = "", x = "") +
  scale_x_discrete(labels = parse(text = c("PM[2.5]", "NO[2]", "Summer~Ozone")))

# # save plot
# ggsave(paste0(dir_results, "parametric_results/all_parametric_results.png"),
#        p,
#        width = 569,
#        height = 368,
#        units = "px")