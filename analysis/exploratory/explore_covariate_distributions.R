rm(list = ls())
gc()

##### 0. Setup #####

library(data.table)
library(fst)
library(purrr)

setDTthreads(threads = 16)

dir_proj <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/"
# dir_data <- paste0(dir_proj, "data/")
ADRD_agg_lagged <- read_fst(paste0(dir_proj, "data/analysis/ADRD_complete_tv.fst"), as.data.table = TRUE)
setnames(ADRD_agg_lagged, old = c("pct_blk", "pct_owner_occ"), new = c("prop_blk", "prop_owner_occ"))
ADRD_agg_lagged[, `:=`(zip = as.factor(zip), year = as.factor(year), cohort = as.factor(cohort), age_grp = as.factor(age_grp), sex = as.factor(sex), race = as.factor(race), dual = as.factor(dual))]

source(paste0(dir_proj, "code/analysis/helper_functions.R"))


##### 1. Explore covariate distributions of aggregated dataset #####

# disaggregate <- function(var){
#   return(unlist(map2(ADRD_agg_lagged[[var]], ADRD_agg_lagged$n_persons, rep))) 
# }

# Not in dataset: AK, HI, PR
table(ADRD_agg_lagged$statecode)

# Check correlations of variables (NOT included in final dataset)
ADRD_agg_lagged[, entry_age := ADRD_age - ADRD_year + ffs_entry_year]
entry_age <- unlist(map2(ADRD_agg_lagged$entry_age, ADRD_agg_lagged$n_persons, rep))
ADRD_age <- unlist(map2(ADRD_agg_lagged$ADRD_age, ADRD_agg_lagged$n_persons, rep))
cor(entry_age, ADRD_age)
set.seed(200)
random_rows <- sample(1:nrow(ADRD_agg_lagged), 10000)
plot(entry_age[random_rows], ADRD_age[random_rows])
# cor(disaggregate("entry_age"), disaggregate("ADRD_age"))
# plot(disaggregate("entry_age"), disaggregate("ADRD_age"))
ADRD_agg_lagged[, entry_age := NULL]

# Check correlations of quantitative confounders included in final dataset
confounder_correlations <- cor(subset(ADRD_agg_lagged_subset, select = zip_quant_var_names))
which((confounder_correlations < 1 & confounder_correlations > 0.5) | confounder_correlations < -0.5, arr.ind = T)

# Explore categorical variables used to stratify
cat("Number of individuals:", sum(ADRD_agg_lagged$n_persons))
cat("Number of male individuals:", sum(ADRD_agg_lagged$n_persons * ADRD_agg_lagged$sexM))
cat("Number of Medicaid-eligible individuals:", sum(ADRD_agg_lagged$n_persons * ADRD_agg_lagged$any_dual))
cat("Number of ADRD hospitalizations:", sum(ADRD_agg_lagged$n_ADRDhosp))
print("Race:")
print(prop.table(table(unlist(map2(ADRD_agg_lagged$race_cat, ADRD_agg_lagged$n_persons, rep)))))
print("Region:")
print(prop.table(table(unlist(map2(ADRD_agg_lagged$region, ADRD_agg_lagged$n_persons, rep)))))
cat("Number of individuals 67-69 yrs old at ADRD hosp/censoring:", sum(ADRD_agg_lagged$n_persons * (ADRD_agg_lagged$ADRD_age >= 67 & ADRD_agg_lagged$ADRD_age <= 69)), "\n")
cat("Number of individuals 70-74 yrs old at ADRD hosp/censoring:", sum(ADRD_agg_lagged$n_persons * (ADRD_agg_lagged$ADRD_age >= 70 & ADRD_agg_lagged$ADRD_age <= 74)), "\n")
cat("Number of individuals 75-79 yrs old at ADRD hosp/censoring:", sum(ADRD_agg_lagged$n_persons * (ADRD_agg_lagged$ADRD_age >= 75 & ADRD_agg_lagged$ADRD_age <= 79)), "\n")
cat("Number of individuals 80-84 yrs old at ADRD hosp/censoring:", sum(ADRD_agg_lagged$n_persons * (ADRD_agg_lagged$ADRD_age >= 80 & ADRD_agg_lagged$ADRD_age <= 84)), "\n")
cat("Number of individuals 85-89 yrs old at ADRD hosp/censoring:", sum(ADRD_agg_lagged$n_persons * (ADRD_agg_lagged$ADRD_age >= 85 & ADRD_agg_lagged$ADRD_age <= 89)), "\n")
cat("Number of individuals 90-94 yrs old at ADRD hosp/censoring:", sum(ADRD_agg_lagged$n_persons * (ADRD_agg_lagged$ADRD_age >= 90 & ADRD_agg_lagged$ADRD_age <= 94)), "\n")
cat("Number of individuals 95+ yrs old at ADRD hosp/censoring:", sum(ADRD_agg_lagged$n_persons * (ADRD_agg_lagged$ADRD_age >= 95)))

# Plot quantitative variables used to stratify
hist(ADRD_agg_lagged$ADRD_year, main="Year of event", xlab="")
hist(ADRD_agg_lagged$ffs_entry_year, main = "FFS entry year", xlab="")
hist(ADRD_agg_lagged$n_years, main = "Number of years in study", xlab="")

# Summarize ZIP-code-level covariates
explore_zip_covs(ADRD_agg_lagged)

hist(ADRD_agg_lagged$mean_bmi, main = "Mean BMI in ZIP code in year before event", xlab="")
hist(ADRD_agg_lagged$smoke_rate, main = "Smoking rate in ZIP code in year before event", xlab="")
hist(ADRD_agg_lagged$medhouseholdincome, main = "Median household income in ZIP code in year before event", xlab="Dollars")
hist(ADRD_agg_lagged$education, main = "Proportion NOT graduated high school in ZIP code in year before event", xlab="")
hist(ADRD_agg_lagged$pct_blk, main = "Proportion Black in ZIP code in year before event", xlab="")
hist(ADRD_agg_lagged$popdensity, main = "Population Density in ZIP code in year before event", xlab="People/Mile^2")
hist(ADRD_agg_lagged$hispanic, main = "Proportion Hispanic in ZIP code in year before event", xlab="")
hist(ADRD_agg_lagged$pct_owner_occ, main = "Proportion of residents who own their home in ZIP code in year before event", xlab="")
hist(ADRD_agg_lagged$poverty, main = "Proportion of residents living below poverty line in year before event", xlab="")
hist(ADRD_agg_lagged$PIR, main = "Median Price-Income Ratio (PIR) in ZIP code in year before event", xlab="")
hist(ADRD_agg_lagged$pr, main = "Mean annual precipitation in year before event", xlab="Inches")

logit <- function(x){
  return(log(x/(1-x)))
}

# Plot transformed covariates, weighting by strata size
hist(sqrt(ADRD_agg_lagged$medhouseholdincome), main = "Sqrt(Median household income in ZIP code in year before event)", xlab="Sqrt(Dollars)")
hist(logit(ADRD_agg_lagged$education), main = "Logit(Proportion NOT graduated high school in ZIP code in year before event)", xlab="") # sqrt transformation would also make distribution symmetric
hist(logit(ADRD_agg_lagged$pct_blk), main = "Logit(Proportion Black in ZIP code in year before event)", xlab="")
hist(logit(ADRD_agg_lagged$hispanic), main = "Logit(Proportion Hispanic in ZIP code in year before event)", xlab="")
hist(log(ADRD_agg_lagged$popdensity), main = "Log(Population Density in ZIP code in year before event)", xlab="Log(People/Mile^2)")
hist(logit(ADRD_agg_lagged$pct_owner_occ), main = "Logit(Proportion of residents who own their home in ZIP code in year before event)", xlab="")
hist(logit(ADRD_agg_lagged$poverty), main = "Logit(Proportion of residents living below poverty line in ZIP code in year before event)", xlab="")
hist(log(ADRD_agg_lagged$PIR), main = "Log(Mean Poverty-Income Ratio (PIR) in ZIP code in year before event)", xlab="") # skew left after log transformation

# Summarize ZIP-level exposures
mean(ADRD_agg_lagged$pm25)
sd(ADRD_agg_lagged$pm25)
mean(ADRD_agg_lagged$no2)
sd(ADRD_agg_lagged$no2)
mean(ADRD_agg_lagged$ozone_summer)
sd(ADRD_agg_lagged$ozone_summer)
cor(ADRD_agg_lagged$pm25, ADRD_agg_lagged$no2)
cor(ADRD_agg_lagged$pm25, ADRD_agg_lagged$ozone_summer)
cor(ADRD_agg_lagged$no2, ADRD_agg_lagged$ozone_summer)

# Save transformed variables in data table
ADRD_agg_lagged[, `:=`(sqrt_medhouseholdincome = sqrt(medhouseholdincome),
                logit_education = logit(education),
                logit_pct_blk = logit(pct_blk),
                logit_hispanic = logit(hispanic),
                logit_pct_owner_occ = logit(pct_owner_occ),
                log_PIR = log(PIR),
                sqrt_no2 = sqrt(no2))]
