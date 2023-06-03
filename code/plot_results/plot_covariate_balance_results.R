library(data.table)
library(ggplot2)

# get directories and classifications of variables
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/git/code/"
source(paste0(dir_code, "constants.R"))

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

# labels
all_cov_bal$Exposure <- factor(all_cov_bal$Exposure, c("PM2.5", "NO2", "Summer Ozone"),
                               c("PM[2.5]", "NO[2]", "Summer~Ozone"))
all_cov_bal$Dataset <- factor(all_cov_bal$Dataset, c("Unadjusted", "Weighted", "Matched"),
                              c("Unadjusted", "GPS Weighting", "GPS Matching"))
covs <- c("education", "hispanic", "mean_bmi", "MW", "NE", "S", "W",
          "PIR", "popdensity", "poverty", "prop_blk", "prop_owner_occ",
          "smoke_rate", "summer_rmax", "summer_tmmx")
cov_names <- c("Education", "% Hispanic", "Body Mass Index", 
               "Region MidWest", "Region NorthEast", "Region South", "Region West",
               "Home Price:Income", "Population Density",
               "Poverty", "% Black", "% Owner Occupied",
               "Smoking Rate", "Summer Humidity", "Summer Temperature")
cov_order <- rev(order(cov_names))
all_cov_bal$Covariate <- 
  factor(all_cov_bal$Covariate, covs[cov_order], cov_names[cov_order])

# plot covariate balance with exposures as facets
ggplot(all_cov_bal, aes(x = AbsoluteCorrelation,
                        y = Covariate,
                        color = Dataset,
                        shape = Dataset,
                        group = Dataset)) +
  facet_grid(~Exposure, labeller = label_parsed) +
  geom_point(size = 4) +
  geom_line(orientation = "y", linewidth = 1) +
  geom_vline(xintercept = 0, linewidth = 2) +
  geom_vline(xintercept = 0.1, linetype = "dashed", 
             linewidth = 1, color = "black") +
  xlab("Absolute Correlation") +
  ylab("") +
  theme_minimal(base_size = 24) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.44)) +
  theme(legend.key.width = unit(1, "cm"), 
        legend.key.height = unit(1, "cm"),
        legend.position = c(.9, .9),
        # axis.text.y = element_text(angle = 30),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"), 
        strip.placement = "outside") +
  labs(color = "", shape = "")
