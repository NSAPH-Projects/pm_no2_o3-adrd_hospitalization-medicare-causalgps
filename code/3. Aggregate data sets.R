# ############################################################################ #
#' Project: Causal ADRD                                     
#' Code: aggregate qid data
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

##### Read denom files ##### 
ffs_entry_exit_adrd <- read_fst(paste0(dir_data, "denom/ffs_entry_exit_adrd.fst"), 
                           as.data.table = TRUE)


##### Create first AD aggregated dataset #####
# no unknown sex or race
AD_agg <- ffs_entry_exit_adrd[sex != 0 & race != 0
                              , .(n_persons = .N, 
                                  n_ADhosp = sum(AD_hosp), 
                                  n_years = first(1 + AD_year - ffs_entry_year)),
                              by = .(ffs_entry_year, AD_zip, AD_age, 
                                     sexM, race_cat, any_dual)]
AD_agg[, sum(n_ADhosp)] # 2929842
write_fst(AD_agg, paste0(dir_data, "aggregated/AD_agg.fst"))

# begin at age 65, no unknown sex or race
AD_agg65 <- ffs_entry_exit_adrd[entry_age == 65 & sex != 0 & race != 0
                              , .(n_persons = .N, 
                                  n_ADhosp = sum(AD_hosp), 
                                  n_years = first(1 + AD_year - ffs_entry_year)),
                              by = .(ffs_entry_year, AD_zip, AD_age, 
                                     sexM, race_cat, any_dual)]
AD_agg65[, sum(n_ADhosp)] # 298947
write_fst(AD_agg65, paste0(dir_data, "aggregated/AD_agg65.fst"))

##### ADRD aggregated dataset #####
# no unknown sex or race
ADRD_agg <- ffs_entry_exit_adrd[sex != 0 & race != 0
                                , .(n_persons = .N, 
                                    n_ADRDhosp = sum(ADRD_hosp), 
                                    n_years = first(1 + ADRD_year - ffs_entry_year)),
                                by = .(ffs_entry_year, ADRD_zip, ADRD_age, 
                                       sexM, race_cat, any_dual)]
ADRD_agg[, sum(n_ADRDhosp)] # 7540591
write_fst(ADRD_agg, paste0(dir_data, "aggregated/ADRD_agg.fst"))


# begin at age 65, no unknown sex or race
ADRD_agg65 <- ffs_entry_exit_adrd[entry_age == 65 & sex != 0 & race != 0
                                , .(n_persons = .N, 
                                    n_ADRDhosp = sum(ADRD_hosp), 
                                    n_years = first(1 + ADRD_year - ffs_entry_year)),
                                by = .(ffs_entry_year, ADRD_zip, ADRD_age, 
                                       sexM, race_cat, any_dual)]
ADRD_agg65[, sum(n_ADRDhosp)] # 968280
write_fst(ADRD_agg65, paste0(dir_data, "aggregated/ADRD_agg65.fst"))
