rm(list = ls())
gc()

##### 0. Setup #####
library(data.table)
library(fst)
library(purrr)
# library(NSAPHutils)

setDTthreads(threads = 16)

dir_data <- "~/nsaph_projects/pm_no2_o3-adrd_hosp-medicare-causalgps/"
ADRD_agg_unlagged <- read_fst(paste0(dir_data, "data/analysis/ADRD_complete.fst"), as.data.table = TRUE)

# Approximate first ADRD shospitalization by requiring no ADRD hosps for 2 years
ADRD_agg <- ADRD_agg_unlagged[ADRD_year - ffs_entry_year >= 2, ]


##### 1. Explore covariate distributions at individual level #####

# Explore categorical variables used to stratify
cat("Number of individuals:", sum(ADRD_agg$n_persons))
cat("Number of male individuals:", sum(ADRD_agg$n_persons * ADRD_agg$sexM))
cat("Number of Medicaid-eligible individuals:", sum(ADRD_agg$n_persons * ADRD_agg$any_dual))
cat("Number of ADRD hospitalizations:", sum(ADRD_agg$n_ADRDhosp))
print("Race:")
print(prop.table(table(unlist(map2(ADRD_agg$race_cat, ADRD_agg$n_persons, rep)))))
print("Region:")
print(prop.table(table(unlist(map2(ADRD_agg$region, ADRD_agg$n_persons, rep)))))

# Plot quantitative variables used to stratify
hist(unlist(map2(ADRD_agg$ADRD_year, ADRD_agg$n_persons, rep)), main="Year of event", xlab="")
hist(unlist(map2(ADRD_agg$ffs_entry_year, ADRD_agg$n_persons, rep)), main = "FFS entry year", xlab="")
hist(unlist(map2(ADRD_agg$n_years, ADRD_agg$n_persons, rep)), main = "Number of years in study", xlab="")

## Plot quantitative covariates, weighting by strata size ##
hist(unlist(map2(ADRD_agg$mean_bmi, ADRD_agg$n_persons, rep)), main = "Mean BMI in ZIP code in year before event", xlab="")
hist(unlist(map2(ADRD_agg$smoke_rate, ADRD_agg$n_persons, rep)), main = "Smoking rate in ZIP code in year before event", xlab="")
hist(unlist(map2(ADRD_agg$medhouseholdincome, ADRD_agg$n_persons, rep)), main = "Median household income in ZIP code in year before event", xlab="Dollars")
hist(unlist(map2(ADRD_agg$education, ADRD_agg$n_persons, rep)), main = "Education rate? in ZIP code in year before event", xlab="")
hist(unlist(map2(ADRD_agg$pct_blk, ADRD_agg$n_persons, rep)), main = "Percent Black in ZIP code in year before event", xlab="")
hist(unlist(map2(ADRD_agg$hispanic, ADRD_agg$n_persons, rep)), main = "Percent Hispanic in ZIP code in year before event", xlab="")
hist(unlist(map2(ADRD_agg$pct_owner_occ, ADRD_agg$n_persons, rep)), main = "Percent Owner Occupancy? in ZIP code in year before event", xlab="")
hist(unlist(map2(ADRD_agg$PIR, ADRD_agg$n_persons, rep)), main = "Mean Poverty-Income Ratio (PIR) in ZIP code in year before event", xlab="")
hist(unlist(map2(ADRD_agg$pr, ADRD_agg$n_persons, rep)), main = "Mean annual precipitation in year before event", xlab="Inches")

logit <- function(x){
  return(log(x/(1-x)))
}

# Plot transformed covariates, weighting by strata size
hist(sqrt(unlist(map2(ADRD_agg$medhouseholdincome, ADRD_agg$n_persons, rep))), main = "Sqrt(Median household income in ZIP code in year before event)", xlab="Sqrt(Dollars)")
hist(logit(unlist(map2(ADRD_agg$education, ADRD_agg$n_persons, rep))), main = "Logit(Education rate? in ZIP code in year before event)", xlab="") # sqrt transformation would also make distribution symmetric
hist(logit(unlist(map2(ADRD_agg$pct_blk, ADRD_agg$n_persons, rep))), main = "Logit(Percent Black in ZIP code in year before event)", xlab="")
hist(logit(unlist(map2(ADRD_agg$hispanic, ADRD_agg$n_persons, rep))), main = "Logit(Percent Hispanic in ZIP code in year before event)", xlab="")
hist(logit(unlist(map2(ADRD_agg$pct_owner_occ, ADRD_agg$n_persons, rep))), main = "Logit(Percent Owner Occupancy? in ZIP code in year before event)", xlab="")
hist(log(unlist(map2(ADRD_agg$PIR, ADRD_agg$n_persons, rep))), main = "Log(Mean Poverty-Income Ratio (PIR) in ZIP code in year before event)", xlab="") # skew left after log transformation

# Plot exposure variables (quantitative), weighted by strata size
hist(unlist(map2(ADRD_agg$pm25, ADRD_agg$n_persons, rep)), main = "Mean annual PM2.5 in ZIP code in year before event", xlab="Micrograms/cubic meter")
hist(unlist(map2(ADRD_agg$no2, ADRD_agg$n_persons, rep)), main = "Mean annual NO2 in ZOP code in year before event", xlab="Micrograms/cubmic meter") # skew right
hist(unlist(map2(ADRD_agg$ozone_summer, ADRD_agg$n_persons, rep)), main = "Mean summer ozone in ZIP code in year before event", xlab="Micrograms/cubic meter")

hist(sqrt(unlist(map2(ADRD_agg$no2, ADRD_agg$n_persons, rep))), main = "Sqrt(Mean annual NO2 in ZOP code in year before event)", xlab="Sqrt(Micrograms/cubmic meter)")

# Save transformed variables in data table
ADRD_agg[, `:=`(sqrt_medhouseholdincome = sqrt(medhouseholdincome),
                logit_education = logit(education),
                logit_pct_blk = logit(pct_blk),
                logit_hispanic = logit(hispanic),
                logit_pct_owner_occ = logit(pct_owner_occ),
                log_PIR = log(PIR),
                sqrt_no2 = sqrt(no2))]

