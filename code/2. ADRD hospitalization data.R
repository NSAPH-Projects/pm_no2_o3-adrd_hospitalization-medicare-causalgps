# ############################################################################ #
#' Project: Causal ADRD                                     
#' Code: combine qid entry/exit with first AD/ADRD hospitalization
#' Inputs: denominator files, QD exposure data, meteorological data, hosp data             
#' Outputs: exposures in grid by qid (rows) and year (cols)                    
#' Author: Daniel Mork                                                         
#' Last updated: Apr 30, 2022                                                  
#' Memory to run: 64 GB
# ############################################################################ #
rm(list = ls())
gc()

##### Setup #####
library(data.table)
library(fst)
library(NSAPHutils)

options(stringsAsFactors = FALSE)
setDTthreads(threads = 16)

dir_data <- "/nfs/home/D/dam9096/shared_space/ci3_analysis/dmork/Data/Causal_ADRD/"
dir_hosp <- "/nfs/home/D/dam9096/shared_space/ci3_analysis/data_ADRDhospitalization/ADRDhospitalization_CCWlist/"
dir_denominator <- "/nfs/home/D/dam9096/shared_space/ci3_health_data/medicare/mortality/1999_2016/wu/cache_data/merged_by_year_v2/"


#### Read AD/ADRD hospitalization data, merge with denominator ####
f <- list.files(dir_hosp, pattern = "\\.fst", full.names = TRUE)
hospvars <- c("QID", "ADATE", "year", "AD_any", "ADRD_any")
ADRDhosp <- rbindlist(lapply(f, read_fst, columns = hospvars, 
                             as.data.table = TRUE))

# First AD hospitalization
setkey(ADRDhosp, QID, ADATE)
firstAnyAD <- ADRDhosp[AD_any == TRUE, 
                       .(AD_date = min(ADATE), AD_year = min(year)), by = QID]
firstAnyAD[, AD_any := TRUE]
setnames(firstAnyAD, "AD_any", "AD_hosp")

# First ADRD hospitalization
firstAnyADRD <- ADRDhosp[ADRD_any == TRUE, 
                         .(ADRD_date = min(ADATE), ADRD_year = min(year)), 
                         by = QID]
firstAnyADRD[, ADRD_any := TRUE]
setnames(firstAnyADRD, "ADRD_any", "ADRD_hosp")
rm(ADRDhosp)

# Combined AD/ADRD first hospitalization
firstHosp <- merge(firstAnyADRD, firstAnyAD, by = "QID", all.x = TRUE)
rm(firstAnyAD, firstAnyADRD)
setkey(firstHosp, QID)


#### Match zipcode to hosp year ####
f <- list.files(dir_denominator, pattern = "\\.fst", full.names = TRUE)
myvars <- c("qid", "year", "zip")
dt <- rbindlist(lapply(f[2:18], read_fst, columns = myvars, as.data.table = TRUE))
setkey(dt, qid, year)
firstHosp <- merge(firstHosp, dt[, .(qid, year, zip)], 
                   by.x = c("QID", "ADRD_year"), by.y = c("qid", "year"), 
                   all.x = TRUE)
firstHosp <- merge(firstHosp, dt[, .(qid, year, zip)], 
                   by.x = c("QID", "AD_year"), by.y = c("qid", "year"), 
                   all.x = TRUE)
setnames(firstHosp, c("zip.x", "zip.y"), c("ADRD_zip", "AD_zip"))

#### Merge into qid entry exit data ####
ffs_entry_exit <- read_fst(paste0(dir_data, "denom/ffs_entry_exit.fst"), 
                           as.data.table = TRUE)
ffs_entry_exit_adrd <- merge(ffs_entry_exit, firstHosp, 
                             by.x = "qid", by.y = "QID", all.x = TRUE)

# data type fixes
ffs_entry_exit_adrd[, c("AD_year", "ADRD_year") := 
                      list(as.integer(AD_year), as.integer(ADRD_year))]

# if no hosp or after FFS censoring, update end year data
ffs_entry_exit_adrd[is.na(AD_hosp) | AD_year > ffs_exit_year, 
                    c("AD_year", "AD_zip", "AD_hosp") :=
                      list(ffs_exit_year, exit_zip, FALSE)]
ffs_entry_exit_adrd[is.na(ADRD_hosp) | ADRD_year > ffs_exit_year,
                    c("ADRD_year", "ADRD_zip", "ADRD_hosp") :=
                      list(ffs_exit_year, exit_zip, FALSE)]

# age at year of event
ffs_entry_exit_adrd[, AD_age := entry_age + AD_year - ffs_entry_year]
ffs_entry_exit_adrd[, ADRD_age := entry_age + ADRD_year - ffs_entry_year] 


##### Categorical Race/Sex Vars #####
ffs_entry_exit_adrd[, sexM := (sex == 1)] # binary sex
race_cats <- c("white", "black", "other", "asian", "hisp", "n_amer_native")
ffs_entry_exit_adrd[race != 0, race_cat := race_cats[race]] # categorical race

# Save
write_fst(ffs_entry_exit_adrd, paste0(dir_data, "denom/ffs_entry_exit_adrd.fst"))
