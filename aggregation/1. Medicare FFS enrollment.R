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
library(NSAPHutils)

options(stringsAsFactors = FALSE)
setDTthreads(threads = 16)

dir_data <- "/n/dominici_nsaph_l3/Lab/projects/pm_no2_o3-adrd_hosp-medicare-causalgps/data/"
dir_denominator <- "/n/dominici_nsaph_l3/Lab/projects/analytic/denom_by_year/"

##### Read denom files ##### 
cat("Reading denominator files...")
f <- list.files(dir_denominator, pattern = "\\.fst", full.names = TRUE)
myvars <- c("qid", "year", "zip", "hmo_mo", "age", "race", "sex", "dual", "dead", "fips")
dt <- rbindlist(lapply(f[2:18], read_fst, columns = myvars, as.data.table = TRUE))
setkey(dt, year, qid) 
dt[,.N]# 651691916 total person years of data
# ! Check how many people removed
dt[year < 2000 | hmo_mo != 0, .N] # 162148450 person years removed

dt <- dt[year >= 2000 & 
           hmo_mo == 0] # retain FFS only (zero HMO months)
setkey(dt, qid, year)
dt <- unique(dt, by = c("qid", "year"))
cat(dt[,.N], "records\n")
# 489543466


##### Enrollment period and continuous enrollment #####
cat("Creating entry/exit data...\n")
# function to get last year of continuous enrollment in Medicare FFS
last_yr_cont <- function(year) 
{ year[first(which(!(first(year):2017 %in% year))) - 1] }
# function to get another variable at last year of enrollment in FFS
last_var_cont <- function(year, var) 
{ var[first(which(!(first(year):2017 %in% year))) - 1] }
# aggregate by qid
ffs_entry_exit <- dt[, .(entry_age = first(age), # get age at entry
                         entry_zip = first(zip), # zip code at enrollment
                         race = first(race), # race should not change
                         sex = first(sex), # sex should not change
                         any_dual = any(dual == "1"), # any dual enrollment, indicator of SE status
                         ffs_entry_year = first(year), # 
                         ffs_exit_year = last_yr_cont(year), # left FFS or died or other censored
                         exit_zip = last_var_cont(year, zip)), 
                     by = qid]
ffs_entry_exit[, .N] # 64404527
setkey(ffs_entry_exit, qid)
write_fst(ffs_entry_exit, paste0(dir_data, "denom/ffs_entry_exit.fst"))
