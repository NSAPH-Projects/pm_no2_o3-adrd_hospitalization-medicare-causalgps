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

dir_data <- "/nfs/home/M/miq336/shared_space/ci3_analysis/dmork/Data/Causal_ADRD/"

##### 0. Load data #####
yr_zip_dat <- read_fst(paste0(dir_data, "denom/year_zip_confounders.fst"), as.data.table = TRUE)
yr_zip_dat[, dat_year := year + 1] # year for merging into dataset (year before end)
yr_zip_dat[, zip := as.integer(zip)]
setkey(yr_zip_dat, zip, dat_year)

expos_dat <- read_fst(paste0(dir_data, "denom/year_zip_exposures.fst"), as.data.table = TRUE)
expos_dat[, dat_year := year + 1] # year for merging into dataset (year before end)
expos_dat[, zip := as.integer(zip)]
setkey(expos_dat, zip, dat_year)

AD_agg <- read_fst(paste0(dir_data, "aggregated/AD_agg.fst"), as.data.table = TRUE)
AD_agg[, AD_year := ffs_entry_year + n_years - 1]
ADRD_agg <- read_fst(paste0(dir_data, "aggregated/ADRD_agg.fst"), as.data.table = TRUE)
ADRD_agg[, ADRD_year := ffs_entry_year + n_years - 1]

#### 1. Merge and clean ####
setkey(AD_agg, AD_zip, AD_year)
AD_agg <- merge(AD_agg, yr_zip_dat, 
                by.x = c("AD_zip", "AD_year"), by.y = c("zip", "dat_year"),
                all.x = TRUE)
AD_agg <- merge(AD_agg, expos_dat, 
                by.x = c("AD_zip", "AD_year"), by.y = c("zip", "dat_year"),
                all.x = TRUE)
setnames(AD_agg, "year.x", "conf_year")

setkey(ADRD_agg, ADRD_zip, ADRD_year)
ADRD_agg <- merge(ADRD_agg, yr_zip_dat, 
                  by.x = c("ADRD_zip", "ADRD_year"), by.y = c("zip", "dat_year"),
                  all.x = TRUE)
ADRD_agg <- merge(ADRD_agg, expos_dat, 
                  by.x = c("ADRD_zip", "ADRD_year"), by.y = c("zip", "dat_year"),
                  all.x = TRUE)
setnames(ADRD_agg, "year.x", "conf_year")

## remove missing zips/statecode
AD_agg <- AD_agg[!is.na(statecode)]
ADRD_agg <- ADRD_agg[!is.na(statecode)]

## remaining complete cases
AD_agg[complete.cases(AD_agg), .N] # 12167271
ADRD_agg[complete.cases(ADRD_agg), .N] # 12094195

## Save complete data
write_fst(AD_agg[complete.cases(AD_agg)], paste0(dir_data, "analysis/AD_complete.fst"))
write_fst(ADRD_agg[complete.cases(ADRD_agg)], paste0(dir_data, "analysis/ADRD_complete.fst"))
