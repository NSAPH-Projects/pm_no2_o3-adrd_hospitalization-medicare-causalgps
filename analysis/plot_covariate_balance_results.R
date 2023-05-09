library(data.table)
library(ggplot2)

# directories for data, code, and results
dir_data <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/data/"
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/code/"
dir_results <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/results/"

# get helpful constants
source(paste0(dir_code, "analysis/constants.R"))

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
all_cov_bal[, Dataset := factor(Dataset, levels = c("Unadjusted", "Matched", "Weighted"))]
all_cov_bal[, Exposure := factor(Exposure, levels = c("pm25", "no2", "ozone_summer"), labels = c("PM2.5", "NO2", "Summer Ozone"))]

# plot covariate balance with exposures as facets
cov_bal_plot <- ggplot(all_cov_bal, aes(x = AbsoluteCorrelation,
                                        y = reorder(Covariate, AbsoluteCorrelation),
                                        color = Dataset,
                                        shape = Dataset,
                                        group = Dataset)) +
  facet_grid(vars(Exposure)) +
  geom_point() +
  geom_line(orientation = "y") +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "black") +
  xlab("Absolute Correlation with Exposure") +
  ylab("Covariate") +
  theme_minimal() +
  theme(legend.position = "bottom")
