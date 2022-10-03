# ############################################################################ #
#' Project: Causal ADRD                                     
#' Code: gather census and BRFSS data for every year/zipcode     
#' Inputs: denominator data (2000-2016), BRFSS data (exclude b/c no zip after 2012)                             
#' Outputs: year/zip confounders   
#' Author: Daniel Mork                                                         
#' Last updated: Mar 09, 2022                                                  
#' Memory to run: <32 GB
# ############################################################################ #
rm(list = ls())
gc()

##### 0. Setup #####
library(data.table)
library(fst)
library(dplyr)
library(NSAPHutils)
library(zipcode)
data(zipcode)
options(stringsAsFactors = FALSE)

setDTthreads(threads = 16)

dir_data <- "/nfs/home/D/dam9096/shared_space/ci3_analysis/dmork/Data/Causal_ADRD/"


##### 1. Create zip-year data from National Causal dataset #####
load(file="/nfs/home/D/dam9096/shared_space/ci3_analysis/National_Causal/National_Causal/2016_temp/aggregate_data.RData")
setDT(aggregate_data)
ad <- aggregate_data[, .(mean_bmi = mean(mean_bmi), 
                         smoke_rate = mean(smoke_rate), 
                         hispanic = mean(hispanic),
                         prop_blk = mean(pct_blk), # between 0 and 1
                         medhouseholdincome = mean(medhouseholdincome), 
                         medianhousevalue = mean(medianhousevalue),
                         PIR = mean(medianhousevalue) / mean(medhouseholdincome),
                         poverty = mean(poverty), 
                         education = mean(education), 
                         popdensity = mean(popdensity),
                         prop_owner_occ = mean(pct_owner_occ), # between 0 and 1
                         summer_tmmx = mean(summer_tmmx), 
                         winter_tmmx = mean(winter_tmmx),
                         summer_rmax = mean(summer_rmax), 
                         winter_rmax = mean(winter_rmax)), 
                     by = .(zip, year)]
ad <- ad[medhouseholdincome > 0]
ad <- merge(ad, zipcode, by = "zip", all.x = TRUE)
setnames(ad, "state", "statecode")


# 4 regions based on census?
NE <- c("NY", "MA", "PA", "RI", "NH", "ME", "VT", "CT", "NJ")  
S <- c("DC", "VA", "NC", "WV", "KY", "SC", "GA", "FL", "AL", "TN", "MS", 
       "AR", "MD", "DE", "OK", "TX", "LA")
MW <- c("OH", "IN", "MI", "IA", "MO", "WI", "MN", "SD", "ND", "IL", "KS", "NE")
W <- c("MT", "CO", "WY", "ID", "UT", "NV", "CA", "OR", "WA", "AZ", "NM")

ad[, region := ifelse(statecode %in% NE, "NE",
                      ifelse(statecode %in% S, "S", 
                             ifelse(statecode %in% MW, "MW",
                                    ifelse(statecode %in% W, "W", "Other"))))]
write_fst(ad, paste0(dir_data, "denom/year_zip_confounders.fst"))



##### 2. Exposure data #####
dir_pm25 <- "/nfs/home/D/dam9096/shared_space/ci3_exposure/pm25/whole_us/annual/zipcode/qd_predictions_ensemble/ywei_aggregation/"
dir_pm25_comp <- "/nfs/home/D/dam9096/shared_space/ci3_exposure/pm25_components/"
dir_no2 <- "/nfs/home/D/dam9096/shared_space/ci3_exposure/no2/whole_us/annual/zipcode/qd_predictions_ensemble/ywei_aggregations/"
dir_ozone <- "/nfs/home/D/dam9096/shared_space/ci3_exposure/ozone/whole_us/seasonal/zipcode/"
dir_temp <- "/nfs/home/D/dam9096/shared_space/ci3_confounders/data_for_analysis/prepped_temperature/annual/"

pm25_data <- fread(paste0(dir_pm25, "all_years.csv"))[, zip := ZIP][, ZIP := NULL][zip >= 501]
no2_data <- fread(paste0(dir_no2, "all_years.csv"))[, zip := ZIP][, ZIP := NULL][zip >= 501]
ozone_data <- fread(paste0(dir_ozone, "ozone_seasonalavg_zipcode.csv"))[, zip := ZIP][, ZIP := NULL][zip >= 501]
temp_rh_data <- fread(paste0(dir_temp, "temperature_annual_zipcode.csv"))[, zip := ZIP][, ZIP := NULL][zip >= 501]
all_exp <- merge(merge(merge(pm25_data, 
                             no2_data, by = c("zip", "year")), 
                       ozone_data, by = c("zip", "year")),
                 temp_rh_data, by = c("zip", "year"))
rm(pm25_data, no2_data, ozone_data, temp_rh_data)
setkey(all_exp, zip, year)

write_fst(all_exp, paste0(dir_data, "denom/year_zip_exposures.fst"))

