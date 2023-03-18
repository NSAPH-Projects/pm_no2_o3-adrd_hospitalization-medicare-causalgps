library(data.table)
library(fst)
library(survival)

# directories for data, code, and results
dir_data <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/data/"
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/code/"
dir_results <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/results/"

# get classifications of variables
source(paste0(dir_code, "analysis/helper_functions.R"))

# get (time-varying) patient data
dt <- read_fst(paste0(dir_data, "denom/complete_ADRD_denom.fst"),
               as.data.table = TRUE)
setkey(dt, zip, cohort, year, age_grp, sex, race, dual)
dt <- dt[ADRD_year >= cohort + 2 & # 2 year 'clean' period
           ADRD_year >= year & # censor observations after ADRD event
           year >= 2001] # begin with year 2001 to pair with 2000 exposures
dt[cohort == 2000, cohort := 2001] # start 2000 cohort in 2001
dt <- subset(dt, select = c("qid",
                            "cohort",
                            "ADRD_year",
                            "last_year_ffs",
                            "ADRD_hosp",
                            "zip",
                            indiv_var_names)) # includes "year"

## Note: No trimming (by exposure or GPS) has been performed at this point

# get survival data
dt_ADRD_year <- dt[year == ADRD_year]
dt_ADRD_year[, n_years := ADRD_year - cohort + 1]

# merge exposures
expos_dat <- read_fst(paste0(dir_data, "denom/year_zip_exposures.fst"),
                      as.data.table = TRUE)
expos_dat[, dat_year := year + 1][, year := NULL] # year for merging into dataset (year before end)
expos_dat[, zip := as.integer(zip)]
setkey(expos_dat, zip, dat_year)
dt_ADRD_year <- merge(dt_ADRD_year, expos_dat,
                       by.x = c("zip", "year"), by.y = c("zip", "dat_year"),
                       all.x = TRUE)

# fit survival curve with sex as predictor
fit <- survfit(Surv(time = n_years, event = ADRD_hosp, type = "right") ~ sex,
               data = dt_ADRD_year)

# save plot
png(paste0(dir_results, "exploratory/km_curve_by_sex.png"))
plot(fit[1],
     main = "Kaplan-Meier Estimate of Survival Function (Associational)",
     xlab = "Number of years observed",
     ylab = "Probability of ADRD hospitalization",
     col = "red") # male in red
lines(fit[2], col = "blue") # female in blue
legend("bottomleft", inset = 0.1, col = c("red", "blue"), legend = c("Male", "Female"), lty = 1)
dev.off()

# fit survival curve with quartilized exposure as predictor, for each exposure
dt_ADRD_year[, pm25_quartiles := cut(pm25, breaks = quantile(pm25, na.rm = T))]
dt_ADRD_year[, no2_quartiles := cut(no2, breaks = quantile(no2, na.rm = T))]
dt_ADRD_year[, ozone_summer_quartiles := cut(ozone_summer, breaks = quantile(ozone_summer, na.rm = T))]

fit_pm25 <- survfit(Surv(time = n_years, event = ADRD_hosp, type = "right") ~ pm25_quartiles,
               data = dt_ADRD_year)
fit_no2 <- survfit(Surv(time = n_years, event = ADRD_hosp, type = "right") ~ no2_quartiles,
                    data = dt_ADRD_year)
fit_ozone_summer <- survfit(Surv(time = n_years, event = ADRD_hosp, type = "right") ~ ozone_summer_quartiles,
                    data = dt_ADRD_year)

# save plots
png(paste0(dir_results, "exploratory/km_curve_by_pm25_quartile.png"))
plot(fit_pm25[1],
     main = "Kaplan-Meier Estimate of Survival Function (Associational)",
     xlab = "Number of years observed",
     ylab = "Probability of ADRD hospitalization",
     col = "red") # lowest level of exposure in red
lines(fit_pm25[2], col = "orange")
lines(fit_pm25[3], col = "green")
lines(fit_pm25[4], col = "blue")
legend <- unique(dt_ADRD_year$pm25_quartiles)[!is.na(unique(dt_ADRD_year$pm25_quartiles))]
legend <- paste(legend, "micrograms/m^3 PM2.5")
legend("bottomleft", inset = 0.1, col = c("red", "orange", "green", "blue"), legend = legend, lty = 1)
dev.off()

png(paste0(dir_results, "exploratory/km_curve_by_no2_quartile.png"))
plot(fit_no2[1],
     main = "Kaplan-Meier Estimate of Survival Function (Associational)",
     xlab = "Number of years observed",
     ylab = "Probability of ADRD hospitalization",
     col = "red") # lowest level of exposure in red
lines(fit_no2[2], col = "orange")
lines(fit_no2[3], col = "green")
lines(fit_no2[4], col = "blue")
legend <- unique(dt_ADRD_year$no2_quartiles)[!is.na(unique(dt_ADRD_year$no2_quartiles))]
legend <- paste(legend, "ppb NO2")
legend("bottomleft", inset = 0.1, col = c("red", "orange", "green", "blue"), legend = legend, lty = 1)
dev.off()

png(paste0(dir_results, "exploratory/km_curve_by_ozone_summer_quartile.png"))
plot(fit_no2[1],
     main = "Kaplan-Meier Estimate of Survival Function (Associational)",
     xlab = "Number of years observed",
     ylab = "Probability of ADRD hospitalization",
     col = "red") # lowest level of exposure in red
lines(fit_no2[2], col = "orange")
lines(fit_no2[3], col = "green")
lines(fit_no2[4], col = "blue")
legend <- unique(dt_ADRD_year$ozone_summer_quartiles)[!is.na(unique(dt_ADRD_year$ozone_summer_quartiles))]
legend <- paste(legend, "ppb Ozone")
legend("bottomleft", inset = 0.1, col = c("red", "orange", "green", "blue"), legend = legend, lty = 1)
dev.off()