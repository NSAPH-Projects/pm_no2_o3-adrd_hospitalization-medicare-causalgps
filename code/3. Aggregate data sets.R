# ############################################################################ #
#' Project: Causal ADRD                                     
#' Code: aggregate qid data
#' Inputs: denominator files, QD exposure data, meteorological data, hosp data             
#' Outputs: exposures in grid by qid (rows) and year (cols)                    
#' Author: Daniel Mork (Michelle Qin edited)                                                     
#' Last updated: June 30, 2022                                                  
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
dir_proj <- "~/nsaph_projects/pm_no2_o3-adrd_hosp-medicare-causalgps/"
dir_data <- paste0(dir_proj, "data/")

##### Read denom files ##### 
ffs_entry_exit_adrd <- read_fst(paste0(dir_data, "denom/ffs_entry_exit_adrd.fst"), 
                           as.data.table = TRUE)


##### Create first AD aggregated dataset #####
# no unknown sex or race
AD_agg <- ffs_entry_exit_adrd[sex != 0 & race != 0
                              , .(n_persons = .N,
                                  n_ADhosp = sum(AD_hosp),
                                  n_years = 1 + AD_year - ffs_entry_year),
                              by = .(ffs_entry_year, AD_zip, AD_age, AD_year,
                                                sexM, race_cat, any_dual)]
AD_agg[, sum(n_ADhosp)] # 2929842
write_fst(AD_agg, paste0(dir_data, "aggregated/AD_agg_corrected.fst"))

# begin at age 65, no unknown sex or race
AD_agg65 <- ffs_entry_exit_adrd[entry_age == 65 & sex != 0 & race != 0
                              , .(n_persons = .N,
                                  n_ADhosp = sum(AD_hosp),
                                  n_years = 1 + AD_year - ffs_entry_year),
                              by = .(ffs_entry_year, AD_zip, AD_age, AD_year,
                                     sexM, race_cat, any_dual)]
AD_agg65[, sum(n_ADhosp)] # 298947
write_fst(AD_agg65, paste0(dir_data, "aggregated/AD_agg65_corrected.fst"))


##### Create first ADRD aggregated dataset #####
# no unknown sex or race
ADRD_agg <- ffs_entry_exit_adrd[sex != 0 & race != 0
                                , .(n_persons = .N, 
                                    n_ADRDhosp = sum(ADRD_hosp), 
                                    n_years = 1 + ADRD_year - ffs_entry_year),
                                by = .(ffs_entry_year, ADRD_zip, ADRD_age, ADRD_year,
                                       sexM, race_cat, any_dual)]
ADRD_agg[, sum(n_ADRDhosp)] # 7540591
write_fst(ADRD_agg, paste0(dir_data, "aggregated/ADRD_agg_corrected.fst"))

# The following are equivalent to ADRD_agg
# ADRD_agg_entry_age <- ffs_entry_exit_adrd[sex != 0 & race != 0
#                                 , .(n_persons = .N, 
#                                     n_ADRDhosp = sum(ADRD_hosp), 
#                                     n_years = 1 + ADRD_year - ffs_entry_year),
#                                 by = .(ffs_entry_year, ADRD_zip, entry_age, ADRD_year,
#                                        sexM, race_cat, any_dual)]
# ADRD_agg_both <- ffs_entry_exit_adrd[sex != 0 & race != 0
#                                 , .(n_persons = .N, 
#                                     n_ADRDhosp = sum(ADRD_hosp), 
#                                     n_years = 1 + ADRD_year - ffs_entry_year),
#                                 by = .(ffs_entry_year, ADRD_zip, ADRD_age, ADRD_year, entry_age,
#                                        sexM, race_cat, any_dual)]

# begin at age 65, no unknown sex or race
ADRD_agg65 <- ffs_entry_exit_adrd[entry_age == 65 & sex != 0 & race != 0
                                , .(n_persons = .N, 
                                    n_ADRDhosp = sum(ADRD_hosp), 
                                    n_years = 1 + ADRD_year - ffs_entry_year),
                                by = .(ffs_entry_year, ADRD_zip, ADRD_age, ADRD_year,
                                       sexM, race_cat, any_dual)]
ADRD_agg65[, sum(n_ADRDhosp)] # 968280
write_fst(ADRD_agg65, paste0(dir_data, "aggregated/ADRD_agg65_corrected.fst"))
