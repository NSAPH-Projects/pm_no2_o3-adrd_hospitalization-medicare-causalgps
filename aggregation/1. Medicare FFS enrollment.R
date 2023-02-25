# ############################################################################ #
#' Project: Causal ADRD                                     
#' Code: determine FFS enrollment period for each individual
#' Inputs: denominator files, QD exposure data, meteorological data, hosp data             
#' Outputs: qid entry and exit from ffs enrollment                 
#' Author: Daniel Mork                                                         
#' Last updated: May 12, 2022                                                  
#' Memory to run: 96 GB
# ############################################################################ #
rm(list = ls())
gc()

##### Setup #####
library(data.table)
library(fst)

options(stringsAsFactors = FALSE)
setDTthreads(threads = 24)
threads_fst(nr_of_threads = 24, reset_after_fork = TRUE)

dir_proj <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/"
dir_data <- "/n/dominici_nsaph_l3/Lab/projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/data/"
dir_denominator <- "/n/dominici_nsaph_l3/Lab/projects/analytic/denom_by_year/"

##### Read denom files ##### 
cat("Reading denominator files...")
f <- list.files(dir_denominator, pattern = "\\.fst", full.names = TRUE)
myvars <- c("qid", "year", "zip", "hmo_mo", "age", "race", "sex", "dual", "dead", "fips")
dt <- rbindlist(lapply(f[2:18], read_fst, columns = myvars, as.data.table = TRUE))
setkey(dt, year, qid) 
dt[,.N] # 651691916 total person years of data
# ! Check how many people removed
dt[year < 2000 | hmo_mo != 0, .N] # 162148450 person years removed

dt <- dt[year >= 2000 & 
           hmo_mo == 0 & # retain FFS only (zero HMO months)
           race != 0 & sex != 0] 
setkey(dt, qid, year)
dt[, hmo_mo := NULL]
dt <- unique(dt, by = c("qid", "year"))
cat(dt[,.N], "records\n") # 489543466


##### Determine last year of continuous enrollment #####
cat("\nDetermine cohort and last year FFS")
# function to get last year of continuous enrollment in Medicare FFS
last_yr_cont <- function(year) 
  year[first(which(!(first(year):2017 %in% year))) - 1]
dt[, cohort := first(year), by = qid]
dt[, last_year_ffs := last_yr_cont(year), by = qid]
# remove records after year without some FFS coverage
dt[year > last_year_ffs, .N] # 12039782 removed
dt <- dt[year <= last_year_ffs]
dt[, .N] # 477503684
# corrected age (in case of inconsistencies)
dt[, age_corrected := first(age) + year - first(year), by = qid]
dt[, age := NULL]


# Merge RTI race code
rti <- rbindlist(lapply(c(2009:2014, 2016), function(y) {
  d <- fread(paste0(dir_proj, "data/auxiliary_medicare_cols/rti_race_", y, ".csv"))
  d[, year := y]
}), fill = TRUE)
rti[, rti_race_cd := as.integer(rti_race_cd)]
rti <- rti[!is.na(rti_race_cd) & rti_race_cd != 0]
dt <- merge(dt, rti, by = c("qid", "year"), all.x = TRUE)
dt[!is.na(rti_race_cd), race := rti_race_cd]
rm(rti)


# Old code summarizing only start/end times of each individual

##### Enrollment period and continuous enrollment #####
# cat("Creating entry/exit data...\n")
# # function to get last year of continuous enrollment in Medicare FFS
# last_yr_cont <- function(year) 
# { year[first(which(!(first(year):2017 %in% year))) - 1] }
# # function to get another variable at last year of enrollment in FFS
# last_var_cont <- function(year, var) 
# { var[first(which(!(first(year):2017 %in% year))) - 1] }
# # aggregate by qid
# ffs_entry_exit <- dt[, .(entry_age = first(age), # get age at entry
#                          entry_zip = first(zip), # zip code at enrollment
#                          race = first(race), # race should not change
#                          sex = first(sex), # sex should not change
#                          any_dual = any(dual == "1"), # any dual enrollment, indicator of SE status
#                          ffs_entry_year = first(year), # 
#                          ffs_exit_year = last_yr_cont(year), # left FFS or died or other censored
#                          exit_zip = last_var_cont(year, zip),
#                          entry_fips = first(fips),
#                          exit_fips = last_var_cont(year, fips)), 
#                      by = qid]
# ffs_entry_exit[, .N] # 64404527
# setkey(ffs_entry_exit, qid)
# write_fst(ffs_entry_exit, paste0(dir_data, "denom/ffs_entry_exit.fst"))
