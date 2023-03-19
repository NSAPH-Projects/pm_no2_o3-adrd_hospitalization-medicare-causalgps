library(data.table)
library(fst)
library(plotrix)
library(ggplot2)

# directories for data, code, and results
dir_data <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/data/"
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/code/"
dir_results <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/results/"

# get classifications of variables
source(paste0(dir_code, "analysis/helper_functions.R"))

# read in full data
ADRD_agg_lagged <- read_fst(paste0(dir_data, "analysis/ADRD_complete_tv.fst"),
                            as.data.table = TRUE)
setnames(ADRD_agg_lagged,
         old = c("pct_blk", "pct_owner_occ"),
         new = c("prop_blk", "prop_owner_occ"))

# get variables at ZIP year level
zip_year_data <- subset(ADRD_agg_lagged,
                        select = c("zip", "year", zip_expos_names, zip_var_names))
zip_year_data <- unique(zip_year_data, by = c("zip", "year"))

# get mean, SD, pairwise correlations of ZIP year exposures
# > mean(zip_year_data$pm25)
# [1] 9.951262
# > sd(zip_year_data$pm25)
# [1] 3.223463
# > mean(zip_year_data$no2)
# [1] 16.91286
# > sd(zip_year_data$no2)
# [1] 9.430883
# > mean(zip_year_data$ozone_summer)
# [1] 45.99236
# > sd(zip_year_data$ozone_summer)
# [1] 7.514933
# > cor(zip_year_data$pm25, zip_year_data$no2)
# [1] 0.3991182
# > cor(zip_year_data$pm25, zip_year_data$ozone_summer)
# [1] 0.2499757
# > cor(zip_year_data$no2, zip_year_data$ozone_summer)
# [1] 0.2503145
# > IQR(zip_year_data$pm25)
# [1] 4.159233
# > IQR(zip_year_data$no2)
# [1] 12.01235
# > IQR(zip_year_data$ozone_summer)
# [1] 9.800937

# plot distribution of exposure across ZIP years
png(paste0(dir_results, "exploratory/zip_pm25_histogram.png"),
    width = 600, height = 400)
hist(zip_year_data$pm25,
     xlab = "PM2.5 (micrograms/m^3)",
     main = paste0("Annual mean PM2.5\nacross ",
                   format(length(unique(zip_year_data$zip)), scientific = F, big.mark = ','),
                   " ZIP codes, 2000-2015")) # recall that patients are matched to exposure of year prior to observation
dev.off()

png(paste0(dir_results, "exploratory/zip_no2_histogram.png"),
    width = 600, height = 400)
hist(zip_year_data$no2,
     xlab = "NO2 (ppb)",
     main = paste0("Annual mean NO2\nacross ",
                   format(length(unique(zip_year_data$zip)), scientific = F, big.mark = ','),
                   " ZIP codes, 2000-2015")) # recall that patients are matched to exposure of year prior to observation
dev.off()

png(paste0(dir_results, "exploratory/zip_ozone_summer_histogram.png"),
    width = 600, height = 400)
hist(zip_year_data$ozone_summer,
     xlab = "Ozone (ppb)",
     main = paste0("Summer mean ozone\nacross ",
                   format(length(unique(zip_year_data$zip)), scientific = F, big.mark = ','),
                   " ZIP codes, 2000-2015")) # recall that patients are matched to exposure of year prior to observation
dev.off()

quantile(zip_year_data$pm25, c(0, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99, 1))
quantile(zip_year_data$no2, c(0, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99, 1))
quantile(zip_year_data$ozone_summer, c(0, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99, 1))

# plot distribution of exposure across patients (can be thought of as patient-weighted ZIP years)
png(paste0(dir_results, "exploratory/patient_pm25_histogram.png"),
    width = 600, height = 400)
weighted.hist(x = ADRD_agg_lagged$pm25,
              w = ADRD_agg_lagged$n_persons,
              breaks = seq(0, max(ADRD_agg_lagged$pm25), by = 2),
              xlab = "PM2.5 (micrograms/m^3)",
              main = paste0("Annual mean PM2.5\nacross ",
                            format(nrow(ADRD_agg_lagged), scientific = F, big.mark = ','),
                            " patients"))
dev.off()

png(paste0(dir_results, "exploratory/patient_no2_histogram.png"),
    width = 600, height = 400)
weighted.hist(x = ADRD_agg_lagged$no2,
              w = ADRD_agg_lagged$n_persons,
              breaks = seq(0, max(ADRD_agg_lagged$no2), by = 5),
              xlab = "NO2 (ppb)",
              main = paste0("Annual mean NO2\nacross ",
                            format(nrow(ADRD_agg_lagged), scientific = F, big.mark = ','),
                            " patients"))
dev.off()

png(paste0(dir_results, "exploratory/patient_ozone_summer_histogram.png"),
    width = 600, height = 400)
weighted.hist(x = ADRD_agg_lagged$ozone_summer,
              w = ADRD_agg_lagged$n_persons,
              breaks = seq(0, max(ADRD_agg_lagged$ozone_summer), by = 5),
              xlab = "Ozone (ppb)",
              main = paste0("Summer mean ozone\nacross ",
                            format(nrow(ADRD_agg_lagged), scientific = F, big.mark = ','),
                            " patients"))
dev.off()

# plot distribution of exposure across ZIP codes
zip_year_summary <- zip_year_data[, .(exposure_year = year - 1, # in our dataset, "year" is year of observation, exposure is taken from previous year
                                      pm25_median = median(pm25),
                                      no2_median = median(no2),
                                      ozone_summer_median = median(ozone_summer),
                                      pm25_0.05 = quantile(pm25, 0.05),
                                      no2_0.05 = quantile(no2, 0.05),
                                      ozone_summer_0.05 = quantile(ozone_summer, 0.05),
                                      pm25_0.95 = quantile(pm25, 0.95),
                                      no2_0.95 = quantile(no2, 0.95),
                                      ozone_summer_0.95 = quantile(ozone_summer, 0.95)),
                                  by = year]

# to do: save as paste0(dir_results, "exploratory/zip_pm25_by_year.png"), etc.
p <- ggplot(zip_year_summary, aes(x = exposure_year,
                                  y = pm25_median)) +
  geom_ribbon(aes(ymin = pm25_0.05,
                  ymax = pm25_0.95),
              fill = "steelblue2") +
  geom_line(color = "red", size = 1) +
  ggtitle("Annual mean PM2.5 across ZIP codes, by year\n(median in red, middle 90% in blue)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("PM2.5 (micrograms/m^3)") +
  xlab("Exposure Year")

p <- ggplot(zip_year_summary, aes(x = exposure_year,
                                  y = no2_median)) +
  geom_ribbon(aes(ymin = no2_0.05,
                  ymax = no2_0.95),
              fill = "steelblue2") +
  geom_line(color = "red", size = 1) +
  ggtitle("Annual mean NO2 across ZIP codes, by year\n(median in red, middle 90% in blue)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("NO2 (ppb)") +
  xlab("Exposure Year")

p <- ggplot(zip_year_summary, aes(x = exposure_year,
                                  y = ozone_summer_median)) +
  geom_ribbon(aes(ymin = ozone_summer_0.05,
                  ymax = ozone_summer_0.95),
              fill = "steelblue2") +
  geom_line(color = "red", size = 1) +
  ggtitle("Summer mean ozone across ZIP codes, by year\n(median in red, middle 90% in blue in blue)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Ozone (ppb)") +
  xlab("Exposure Year")