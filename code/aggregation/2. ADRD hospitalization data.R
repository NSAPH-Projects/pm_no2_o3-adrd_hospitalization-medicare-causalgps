# ############################################################################ #
#' Project: Causal ADRD                                     
#' Code: combine qid entry/exit with first AD/ADRD hospitalization
#' Inputs: denominator files, QD exposure data, meteorological data, hosp data             
#' Outputs: exposures in grid by qid (rows) and year (cols)                    
#' Author: Daniel Mork (Michelle Qin edited)                                    
#' Last updated: Oct 18, 2022                                                  
#' Memory to run: 96 GB
# ############################################################################ #
# rm(list = ls())
# gc()
# 
# ##### Setup #####
library(data.table)
library(fst)
# library(NSAPHutils)

# options(stringsAsFactors = FALSE)
# setDTthreads(threads = 16)

# get directories and classifications of variables
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/git/code/"
source(paste0(dir_code, "constants.R"))

source(paste0(dir_code, "aggregation/1. Medicare FFS enrollment.R")) # run previous script
# dir_hosp <- "/n/dominici_nsaph_l3/Lab/projects/analytic/adrd_hospitalization/" # included in code/constants.R


#### Read AD/ADRD hospitalization data, merge with denominator ####
cat("\nRead AD/ADRD hospitalization data")
f <- list.files(dir_hosp, pattern = "\\.fst", full.names = TRUE)
hospvars <- c("QID", "ADATE", "year", "AD_any", "ADRD_any")
ADRDhosp <- rbindlist(lapply(f, read_fst, columns = hospvars, 
                             as.data.table = TRUE))

# First AD hospitalization
setkey(ADRDhosp, QID, ADATE)
firstAnyAD <- ADRDhosp[AD_any == TRUE, .(AD_year = min(year)), by = QID]
firstAnyAD[, AD_any := TRUE]
setnames(firstAnyAD, "AD_any", "AD_hosp")

# First ADRD hospitalization
firstAnyADRD <- ADRDhosp[ADRD_any == TRUE, .(ADRD_year = min(year)), 
                         by = QID]
firstAnyADRD[, ADRD_any := TRUE]
setnames(firstAnyADRD, "ADRD_any", "ADRD_hosp")
rm(ADRDhosp)

# Match firstHosp to qid/year in dt
cat("\nMerge hospitalization data")
dt <- merge(dt, firstAnyAD, by.x = "qid", by.y = "QID", all.x = TRUE)
dt <- merge(dt, firstAnyADRD, by.x = "qid", by.y = "QID", all.x = TRUE)
dt[is.na(AD_hosp), AD_hosp := FALSE]
dt[is.na(ADRD_hosp), ADRD_hosp := FALSE]
dt[(AD_hosp), AD_hosp := (year == AD_year)]
dt[(ADRD_hosp), ADRD_hosp := (year == ADRD_year)]
dt[is.na(AD_year), AD_year := last_year_ffs]
dt[is.na(ADRD_year), ADRD_year := last_year_ffs]
rm(firstAnyAD, firstAnyADRD)


dt[, age_grp := cut(age_corrected, breaks = c(64, 69, 74, 79, 84, 89, 94, Inf))]

write_fst(dt, paste0(dir_data, "denom/complete_ADRD_denom.fst"))
# Old code


# # Combined AD/ADRD first hospitalization
# firstHosp <- merge(firstAnyADRD, firstAnyAD, by = "QID", all.x = TRUE)
# rm(firstAnyAD, firstAnyADRD)
# setkey(firstHosp, QID)
# 
# #### Match zipcode to hosp year ####
# f <- list.files(dir_denominator, pattern = "\\.fst", full.names = TRUE)
# myvars <- c("qid", "year", "zip")
# dt <- rbindlist(lapply(f[2:18], read_fst, columns = myvars, as.data.table = TRUE))
# setkey(dt, qid, year)
# firstHosp <- merge(firstHosp, dt[, .(qid, year, zip)], 
#                    by.x = c("QID", "ADRD_year"), by.y = c("qid", "year"), 
#                    all.x = TRUE)
# firstHosp <- merge(firstHosp, dt[, .(qid, year, zip)], 
#                    by.x = c("QID", "AD_year"), by.y = c("qid", "year"), 
#                    all.x = TRUE)
# setnames(firstHosp, c("zip.x", "zip.y"), c("ADRD_zip", "AD_zip"))
# 
# #### Merge into qid entry exit data ####
# ffs_entry_exit <- read_fst(paste0(dir_data, "denom/ffs_entry_exit.fst"), 
#                            as.data.table = TRUE)
# ffs_entry_exit_adrd <- merge(ffs_entry_exit, firstHosp, 
#                              by.x = "qid", by.y = "QID", all.x = TRUE)
# 
# # data type fixes
# ffs_entry_exit_adrd[, c("AD_year", "ADRD_year") := 
#                       list(as.integer(AD_year), as.integer(ADRD_year))]
# 
# # if no hosp or after FFS censoring, update end year data
# ffs_entry_exit_adrd[is.na(AD_hosp) | AD_year > ffs_exit_year, 
#                     c("AD_year", "AD_zip", "AD_hosp") :=
#                       list(ffs_exit_year, exit_zip, FALSE)]
# ffs_entry_exit_adrd[is.na(ADRD_hosp) | ADRD_year > ffs_exit_year,
#                     c("ADRD_year", "ADRD_zip", "ADRD_hosp") :=
#                       list(ffs_exit_year, exit_zip, FALSE)]
# 
# # age at year of event
# ffs_entry_exit_adrd[, AD_age := entry_age + AD_year - ffs_entry_year]
# ffs_entry_exit_adrd[, ADRD_age := entry_age + ADRD_year - ffs_entry_year] 
# 
# 
# ##### Categorical Race/Sex Vars #####
# ffs_entry_exit_adrd[, sexM := (sex == 1)] # binary sex
# race_cats <- c("white", "black", "other", "asian", "hisp", "n_amer_native")
# ffs_entry_exit_adrd[race != 0, race_cat := race_cats[race]] # categorical race
# 
# # Save
# write_fst(ffs_entry_exit_adrd, paste0(dir_data, "denom/ffs_entry_exit_adrd.fst"))
