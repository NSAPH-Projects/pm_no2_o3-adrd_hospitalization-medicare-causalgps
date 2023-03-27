library(data.table)
library(ggplot2)

# directories for data, code, and results
dir_data <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/data/"
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/code/"
dir_results <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/results/"

# get helpful constants
source(paste0(dir_code, "analysis/constants.R"))

# # set up data.table to contain all analyses' covariate balance results
# cov_bal <- expand.grid(Exposure = zip_expos_names,
#                        Dataset = c("Unadjusted", "GPS Weighting", "GPS Matching"),
#                        Covariate = zip_var_names,
#                        AbsoluteCorrelation = -100.1) # placeholder value; a double

# get covariate balance results from tuned modeling
pm25_weighting <- fread(paste0(dir_results, "covariate_balance/pm25/weighting/attempt3/best_cov_bal.csv"))
no2_weighting <- fread(paste0(dir_results, "covariate_balance/no2/weighting/attempt20/best_cov_bal.csv"))
ozone_summer_weighting <- fread(paste0(dir_results, "covariate_balance/ozone_summer/weighting/attempt27/best_cov_bal.csv"))
pm25_matching <- fread(paste0(dir_results, "covariate_balance/pm25/matching/match_zips_gps_untrimmed_caliper2_attempt23/best_cov_bal.csv"))
no2_matching <- fread(paste0(dir_results, "covariate_balance/no2/matching/match_zips_gps_untrimmed_caliper3.5_attempt8/best_cov_bal.csv"))
ozone_summer_matching <- fread(paste0(dir_results, "covariate_balance/ozone_summer/matching/match_zips_gps_untrimmed_caliper4.5_attempt20/best_cov_bal.csv"))

# label the data.tables, in preparation to merge
pm25_weighting$Exposure <- "pm25"
no2_weighting$Exposure <- "no2"
ozone_summer_weighting$Exposure <- "ozone_summer"
pm25_matching$Exposure <- "pm25"
no2_matching$Exposure <- "no2"
ozone_summer_matching$Exposure <- "ozone_summer"

# combine the data.tables
all_cov_bal <- rbindlist(list(pm25_weighting, no2_weighting, ozone_summer_weighting,
                              pm25_matching, no2_matching, ozone_summer_matching))
all_cov_bal <- all_cov_bal[Dataset != "Unmatched"] # remove duplicate rows: "Unmatched" and "Unweighted" are the same
all_cov_bal[, Dataset := ifelse(Dataset == "Unweighted", "Unadjusted", Dataset)]
# all_cov_bal <- all_cov_bal[, .(Exposure, Dataset, Covariate, AbsoluteCorrelation)]
# all_cov_bal <- dcast(all_cov_bal, Exposure + Covariate ~ Dataset, value.var = "AbsoluteCorrelation")

# plot covariate balance for each exposure
pm25_cov_bal_plot <- ggplot(all_cov_bal[Exposure == "pm25"], aes(x = AbsoluteCorrelation,
                                                                 y = reorder(Covariate, AbsoluteCorrelation),
                                                                 color = Dataset,
                                                                 group = Dataset)) +
  geom_point() +
  geom_line(orientation = "y") +
  xlab("Absolute Correlation with PM2.5") +
  ylab("Covariate")

no2_cov_bal_plot <- ggplot(all_cov_bal[Exposure == "no2"], aes(x = AbsoluteCorrelation,
                                                               y = reorder(Covariate, AbsoluteCorrelation),
                                                               color = Dataset,
                                                               group = Dataset)) +
  geom_point() +
  geom_line(orientation = "y") +
  xlab("Absolute Correlation with NO2") +
  ylab("Covariate")

ozone_summer_cov_bal_plot <- ggplot(all_cov_bal[Exposure == "ozone_summer"], aes(x = AbsoluteCorrelation,
                                                                                 y = reorder(Covariate, AbsoluteCorrelation),
                                                                                 color = Dataset,
                                                                                 group = Dataset)) +
  geom_point() +
  geom_line(orientation = "y") +
  xlab("Absolute Correlation with Summer Ozone") +
  ylab("Covariate")
