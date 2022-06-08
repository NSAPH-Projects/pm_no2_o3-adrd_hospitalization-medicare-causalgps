rm(list = ls())
gc()

##### 0. Setup #####
library(data.table)
install.packages("fst")
library(fst)
library(purrr)
# library(NSAPHutils)
install.packages("devtools")
library("devtools")
install_github("fasrc/CausalGPS")
library("CausalGPS")

setDTthreads(threads = 16)

dir_data <- "/n/home_fasse/mqin/nsaph_projects/pm_no2_o3-adrd_hosp-medicare-causalgps/"
ADRD_agg_unlagged <- read_fst(paste0(dir_data, "data/analysis/ADRD_complete.fst"), as.data.table = TRUE)

# > names(ADRD_agg)
# [1] "AD_zip"           "ADRD_year"          "ffs_entry_year"    
# [4] "ADRD_age"           "entry_sex"          "entry_race"        
# [7] "exit_dual"            "n_persons"          "n_ADRD_hosp"        
# [10] "n_years"            "conf_year"          "mean_bmi"          
# [13] "smoke_rate"         "hispanic"           "pct_blk"           
# [16] "medhouseholdincome" "medianhousevalue"   "poverty"           
# [19] "education"          "popdensity"         "pct_owner_occ"     
# [22] "summer_tmmx"        "winter_tmmx"        "summer_rmax"       
# [25] "winter_rmax"        "city"               "statecode"         
# [28] "latitude"           "longitude"          "pm25"              
# [31] "no2"                "ozone_summer"       "ozone_winter"      
# [34] "ozone_fall"         "ozone_spring"       "tmmx"              
# [37] "rmax"               "pr"  

# Approximate first ADRD shospitalization by requiring no ADRD hosps for 2 years
ADRD_agg <- ADRD_agg_unlagged[ADRD_year - ffs_entry_year >= 2, ]


##### 1. Explore covariate distributions at individual level #####

# Categorical variables
cat("Number of individuals:", sum(ADRD_agg$n_persons))
cat("Number of male individuals:", sum(ADRD_agg$n_persons * ADRD_agg$sexM))
cat("Number of Medicaid-eligible individuals:", sum(ADRD_agg$n_persons * ADRD_agg$any_dual))
cat("Number of ADRD hospitalizations:", sum(ADRD_agg$n_ADRDhosp))
print("Race:")
print(table(unlist(map2(ADRD_agg$race_cat, ADRD_agg$n_persons, rep)))) # To Do: get names to show up in table
print("Region:")
print(table(unlist(map2(ADRD_agg$region, ADRD_agg$n_persons, rep))))

# Quantitative variables
hist(unlist(map2(ADRD_agg$ADRD_year, ADRD_agg$n_persons, rep)), main="Year of event", xlab="")
hist(unlist(map2(ADRD_agg$ffs_entry_year, ADRD_agg$n_persons, rep)), main = "FFS entry year", xlab="")
hist(unlist(map2(ADRD_agg$n_years, ADRD_agg$n_persons, rep)), main = "Number of years in study", xlab="")
hist(unlist(map2(ADRD_agg$mean_bmi, ADRD_agg$n_persons, rep)), main = "Mean BMI in ZIP code", xlab="")
hist(unlist(map2(ADRD_agg$smoke_rate, ADRD_agg$n_persons, rep)), main = "Smoking rate in ZIP code", xlab="")
hist(unlist(map2(ADRD_agg$medhouseholdincome, ADRD_agg$n_persons, rep)), main = "Median household income in ZIP code", xlab="Dollars")
hist(unlist(map2(ADRD_agg$education, ADRD_agg$n_persons, rep)), main = "Education rate? in ZIP code", xlab="")

# Exposure variables (quantitative)
hist(unlist(map2(ADRD_agg$pm25, ADRD_agg$n_persons, rep)), main = "Mean annual PM2.5", xlab="micrograms/cubic meter")
hist(unlist(map2(ADRD_agg$no2, ADRD_agg$n_persons, rep)), main = "Mean annual NO2", xlab="micrograms/cubmic meter")
hist(unlist(map2(ADRD_agg$ozone_summer, ADRD_agg$n_persons, rep)), main = "Mean summer ozone", xlab="micrograms/cubic meter")
hist(unlist(map2(ADRD_agg$pr, ADRD_agg$n_persons, rep)), main = "Mean annual precipitation", xlab="units")

