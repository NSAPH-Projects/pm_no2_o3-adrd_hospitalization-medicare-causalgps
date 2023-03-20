library(data.table)
library(ggplot2)

# directories for data, code, and results
dir_data <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/data/"
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/code/"
dir_results <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/results/"

# get helpful constants
source(paste0(dir_code, "analysis/constants.R"))

# # get results
# parametric_results_table <- fread(paste0(dir_results, "parametric_results/parametric_results_table.csv"))
# 
# # get 95% confidence intervals for regression coefficient
# parametric_results_table[, `:=`(coef_95CI_lower = coefficient - 1.96 * se_bootstrap,
#                                 coef_95CI_upper = coefficient + 1.96 * se_bootstrap)]
# 
# # convert coefficient to hazard ratio
# parametric_results_table$HR <- 0.1 # placeholder value; a double
# parametric_results_table[exposure == "pm25", HR := exp(4 * coefficient)]
# parametric_results_table[exposure == "no2", HR := exp(10 * coefficient)]
# parametric_results_table[exposure == "ozone_summer", HR := exp(10 * coefficient)]
# 
# # get 95% confidence intervals for HRs
# parametric_results_table[exposure == "pm25", `:=`(HR_95CI_lower = exp(4 * coef_95CI_lower),
#                                                   HR_95CI_upper = exp(4 * coef_95CI_upper))]
# parametric_results_table[exposure == "no2", `:=`(HR_95CI_lower = exp(10 * coef_95CI_lower),
#                                                  HR_95CI_upper = exp(10 * coef_95CI_upper))]
# parametric_results_table[exposure == "ozone_summer", `:=`(HR_95CI_lower = exp(10 * coef_95CI_lower),
#                                                           HR_95CI_upper = exp(10 * coef_95CI_upper))]

# for now, hard code results from google sheet
parametric_results_table <- expand.grid(Exposure = toupper(zip_expos_names), # note: important to NOT call this column "exposure_name", or else won't be able to match with results later
                                        Method = c("Associational", "GPS Weighting", "GPS Matching"),
                                        HR = -1000.1,
                                        HR_95CI_lower = -1000.1,
                                        HR_95CI_upper = -1000.1) # these are placeholder numbers (important that they are doubles)
parametric_results_table <- setorder(parametric_results_table, Exposure)
parametric_results_table$HR <- c(1.14604857,
                                 1.103757204,
                                 1.106011167,
                                 1.118960355,
                                 1.04780762,
                                 1.002202422,
                                 1.057883274,
                                 1.045808217,
                                 1.021528464)
parametric_results_table$HR_95CI_lower <- c(1.143767749,
                                            1.092011656,
                                            1.078081138,
                                            1.117171727,
                                            1.042926089,
                                            0.989966997,
                                            1.05581206,
                                            1.037015391,
                                            1.008600057)
parametric_results_table$HR_95CI_upper <- c(1.14833394,
                                            1.115629086,
                                            1.134664784,
                                            1.120751847,
                                            1.052711999,
                                            1.014589069,
                                            1.05995855,
                                            1.054675598,
                                            1.03462259)

# plot results
p <- ggplot(parametric_results_table, aes(x = Method,
                                          y = HR,
                                          color = Exposure)) +
  geom_point(position=position_dodge(0.75)) +
  geom_errorbar(aes(ymin = HR_95CI_lower,
                    ymax = HR_95CI_upper),
                position=position_dodge(0.75), width = 0.5)
p

# # save plot
# ggsave(paste0(dir_results, "parametric_results/all_parametric_results.png"),
#        p,
#        width = 569,
#        height = 368,
#        units = "px")
