rm(list = ls())
gc()

##### 0. Setup #####
library(data.table)
library(fst)
library(dplyr)
options(stringsAsFactors = FALSE)

setDTthreads(threads = 16)
dir_proj <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/"
# dir_data <- paste0(dir_proj, "data/")

##### 0. Load data #####
yr_zip_dat <- read_fst(paste0(dir_proj, "data/denom/year_zip_confounders.fst"), as.data.table = TRUE)
yr_zip_dat[, dat_year := year + 1][, year := NULL] # year for merging into dataset (year before end)
yr_zip_dat[, zip := as.integer(zip)]
setkey(yr_zip_dat, zip, dat_year)

expos_dat <- read_fst(paste0(dir_proj, "data/denom/year_zip_exposures.fst"),
                      as.data.table = TRUE)
expos_dat[, dat_year := year + 1][, year := NULL] # year for merging into dataset (year before end)
expos_dat[, zip := as.integer(zip)]
setkey(expos_dat, zip, dat_year)

AD_agg <- read_fst(paste0(dir_proj, "data/aggregated/AD_agg_tv.fst"),
                   as.data.table = TRUE)
ADRD_agg <- read_fst(paste0(dir_proj, "data/aggregated/ADRD_agg_tv.fst"),
                     as.data.table = TRUE)

#### 1. Merge and clean ####
setkey(AD_agg, zip, year)
AD_agg <- merge(AD_agg, yr_zip_dat, 
                by.x = c("zip", "year"), by.y = c("zip", "dat_year"),
                all.x = TRUE)
AD_agg <- merge(AD_agg, expos_dat, 
                by.x = c("zip", "year"), by.y = c("zip", "dat_year"),
                all.x = TRUE)

setkey(ADRD_agg, zip, year)
ADRD_agg <- merge(ADRD_agg, yr_zip_dat, 
                  by.x = c("zip", "year"), by.y = c("zip", "dat_year"),
                  all.x = TRUE)
ADRD_agg <- merge(ADRD_agg, expos_dat, 
                  by.x = c("zip", "year"), by.y = c("zip", "dat_year"),
                  all.x = TRUE)

## remove missing zips/statecode
AD_agg <- AD_agg[!is.na(statecode)]
ADRD_agg <- ADRD_agg[!is.na(statecode)]

## remaining complete cases
AD_agg[complete.cases(AD_agg), .N] # 35490132
ADRD_agg[complete.cases(ADRD_agg), .N] # 34763397

## Save complete data
write_fst(AD_agg[complete.cases(AD_agg)], paste0(dir_proj, "data/analysis/AD_complete_tv.fst"))
write_fst(ADRD_agg[complete.cases(ADRD_agg)], paste0(dir_proj, "data/analysis/ADRD_complete_tv.fst"))
