###############################################################################
# Project: Causal effect of air pollution on first ADRD hospitalization       #
# Code: merge exposure data with hospitalization and denominator data         #
# Input: "ADRD_aggregated_entry_to_followup.fst"                              #
# Input: "zip_year_census_data.fst"                                           #
# Input: pm25, no2, and ozone data                                            #                             
# Output: "ADRD_aggregated_entry_to_followup_merged.fst"                      #
# Author: Daniel Mork                                                         #
# Date: 2021-12-29                                                            #
###############################################################################
rm(list = ls())
gc()
############################# 0. Setup ########################################
library(NSAPHutils)
library(data.table)
library(fst)
setDTthreads(threads = 16)

dir_pm25 <- "/nfs/home/D/dam9096/shared_space/ci3_exposure/pm25/whole_us/annual/zipcode/qd_predictions_ensemble/ywei_aggregation/"
dir_no2 <- "/nfs/home/D/dam9096/shared_space/ci3_exposure/no2/whole_us/annual/zipcode/qd_predictions_ensemble/ywei_aggregations/"
dir_ozone <- "/nfs/home/D/dam9096/shared_space/ci3_exposure/ozone/whole_us/annual/zipcode/requaia_predictions/ywei_aggregation/"
dir_data <- "/nfs/home/D/dam9096/shared_space/ci3_analysis/dmork/Data/Causal_ADRD/"

############################# 1. merge ########################################
pm25_data <- fread(paste0(dir_pm25, "all_years.csv"))
pm25_data[, zip2 := int_to_zip_str(ZIP)][, ZIP := NULL]
no2_data <- fread(paste0(dir_no2, "all_years.csv"))
no2_data[, zip2 := int_to_zip_str(ZIP)][, ZIP := NULL]
ozone_data <- fread(paste0(dir_ozone, "all_years.csv"))
ozone_data[, zip2 := int_to_zip_str(ZIP)][, ZIP := NULL]

exposure <- merge(pm25_data, no2_data, by = c("zip2", "year"))
exposure <- merge(exposure, ozone_data, by = c("zip2", "year"))
exposure[, nox := (1.07*no2 + 2.075*ozone)/3.14]
head(exposure)

rm(pm25_data)
rm(no2_data)
rm(ozone_data)
gc()

###################### 2. Merge with zip-year census and ADRD data #############
ADRD_aggregated <- read_fst(paste0(dir_data, "ADRD_aggregated_entry_to_followup.fst"), 
                            as.data.table = T)
zipYrDat <- read_fst(paste0(dir_data, "zip_year_census_data.fst"), as.data.table = T)
ADRD_aggregated <- merge(ADRD_aggregated, zipYrDat[, -c("N")], 
                         by.x = c("followupZip", "followupYear", "sex", "race", "dual"),
                         by.y = c("zip", "year", "sex", "race", "dual"), 
                         all.x = TRUE, all.y = FALSE, allow.cartesian = FALSE)
rm(zipYrDat)
ADRD_aggregated <- merge(ADRD_aggregated, exposure, 
                         by.x = c("followupZip", "followupYear"), 
                         by.y = c("zip2", "year"), 
                         all.x = TRUE, all.y = FALSE)
rm(exposure)

####################### 3. Merge state and region ##############################
library(zipcode)
NE <- c("NY", "MA", "PA", "RI", "NH", "ME", "VT", "CT", "NJ")  
S <- c("DC", "VA", "NC", "WV", "KY", "SC", "GA", "FL", "AL", "TN", "MS", 
       "AR", "MD", "DE", "OK", "TX", "LA")
MW <- c("OH", "IN", "MI", "IA", "MO", "WI", "MN", "SD", "ND", "IL", "KS", "NE")
W <- c("MT", "CO", "WY", "ID", "UT", "NV", "CA", "OR", "WA", "AZ", "NM")
data(zipcode)
ADRD_aggregated <- merge(ADRD_aggregated, zipcode, 
                         by.x = "followupZip", by.y = "zip", 
                         all.x = TRUE, all.y = FALSE)
rm(zipcode)
ADRD_aggregated[, region := ifelse(state %in% NE, "NE", 
                                   ifelse(state %in% S, "S", 
                                          ifelse(state %in% MW, "MW", "W")))]
rm(NE, S, MW, W)
write_fst(ADRD_aggregated, paste0(dir_data, "ADRD_aggregated_entry_to_followup_merged.fst"))
