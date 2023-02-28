# ############################################################################ #
#' Project: Causal ADRD                                     
#' Code: aggregate qid data
#' Inputs: denominator files, QD exposure data, meteorological data, hosp data             
#' Outputs: exposures in grid by qid (rows) and year (cols)                    
#' Author: Daniel Mork (Michelle Qin edited)                                                     
#' Last updated: Oct 18, 2022                                                  
#' Memory to run: 128 GB
# ############################################################################ #
# rm(list = ls())
# gc()
# 
# ##### Setup #####
# library(data.table)
# library(fst)
# library(NSAPHutils)
# options(stringsAsFactors = FALSE)
# 
# setDTthreads(threads = 16)
dir_proj <- "~/nsaph_projects/pm_no2_o3-adrd_hosp-medicare-causalgps/"
source(paste0(dir_proj, "code/aggregation/2. ADRD hospitalization data.R"))

##### Read denom files ##### 
dt <- read_fst(paste0(dir_proj, "data/denom/complete_ADRD_denom.fst"), as.data.table = TRUE)
setkey(dt, zip, cohort, year, age_grp, sex, race, dual)


# Create ADRD aggregation:
#   Limit: 2 year 'clean' period, stop following after event year, 
#         begin year 2001 to pair with 2000 exposures
#   Select: # persons, # hospitalizations
#   By: zip, cohort, year, age_grp, sex, race, dual
dt_ADRD <- dt[ADRD_year >= cohort + 2 & ADRD_year >= year & year >= 2001]
dt_ADRD[cohort == 2000, cohort := 2001] # start 2000 cohort in 2001
ADRD_agg <- dt_ADRD[, .(n_persons = .N, n_hosp = sum(ADRD_hosp)),
                    by = .(zip, cohort, year, age_grp, sex, race, dual)]
ADRD_agg[, n_years := 1 + year - cohort]
ADRD_agg[, .(.N, sum(n_hosp))] # 41948558 units 5935558 events
write_fst(ADRD_agg, paste0(dir_proj, "data/aggregated/ADRD_agg_tv.fst")) 
rm(dt_ADRD)

# Create AD aggregation:
#   Limit: 2 year 'clean' period, stop following after event year, 
#         begin year 2001 to pair with 2000 exposures
#   Select: # persons, # hospitalizations
#   By: zip, cohort, year, age_grp, sex, race, dual
dt_AD <- dt[AD_year >= cohort + 2 & AD_year >= year & year >= 2001]
dt_AD[cohort == 2000, cohort := 2001] # start 2000 cohort in 2001
AD_agg <- dt_AD[, .(n_persons = .N, n_hosp = sum(AD_hosp)),
                by = .(zip, cohort, year, age_grp, sex, race, dual)]
AD_agg[, n_years := 1 + year - cohort]
AD_agg[, .(.N, sum(n_hosp))] # 42797243 units 2325106 events
write_fst(AD_agg, paste0(dir_proj, "data/aggregated/AD_agg_tv.fst")) 


# ffs_entry_exit_adrd <- read_fst(paste0(dir_proj, "data/denom/ffs_entry_exit_adrd.fst"), 
#                            as.data.table = TRUE)

##### Bin ages by (mostly) 5-year categories #####
# max_AD_age <- max(ffs_entry_exit_adrd$AD_age)
# max_ADRD_age <- max(ffs_entry_exit_adrd$ADRD_age)
# ffs_entry_exit_adrd[, AD_age_binned := cut(AD_age, breaks = c(64, 69, 74, 79, 84, 89, 94, max_AD_age))]
# ffs_entry_exit_adrd[, ADRD_age_binned := cut(ADRD_age, breaks = c(64, 69, 74, 79, 84, 89, 94, max_AD_age))]


##### Create first AD aggregated dataset, binning ages #####
# AD_agg <- dt[sex != 0 & race != 0 & year <= AD_year,
#              .(n_persons = .N, n_hosp = sum(AD_hosp), n_years = 1 + year - cohort),
#              by = .(zip, cohort, year, age_grp, sex, race, dual)]
# AD_agg[, .(.N, sum(n_hosp), sum(n_persons), sum(n_persons * n_years))] 
# # 50,473,844 units, 2907440 events, 
# write_fst(AD_agg, paste0(dir_proj, "data/aggregated/AD_agg_tv.fst")) # AD_agg_corrected.fst was pre-age-binning

# no unknown sex or race
# AD_agg <- ffs_entry_exit_adrd[sex != 0 & race != 0
#                               , .(n_persons = .N,
#                                   n_ADhosp = sum(AD_hosp),
#                                   n_years = 1 + AD_year - ffs_entry_year),
#                               by = .(ffs_entry_year, AD_zip, AD_age_binned, AD_year,
#                                                 sexM, race_cat, any_dual)]
# AD_agg[, sum(n_ADhosp)] # 2929842
# write_fst(AD_agg, paste0(dir_proj, "data/aggregated/AD_agg_age_binned.fst")) # AD_agg_corrected.fst was pre-age-binning

# begin at age 65, no unknown sex or race (pre-age-binning)
# AD_agg65 <- ffs_entry_exit_adrd[entry_age == 65 & sex != 0 & race != 0
#                               , .(n_persons = .N,
#                                   n_ADhosp = sum(AD_hosp),
#                                   n_years = 1 + AD_year - ffs_entry_year),
#                               by = .(ffs_entry_year, AD_zip, AD_age, AD_year,
#                                      sexM, race_cat, any_dual)]
# AD_agg65[, sum(n_ADhosp)] # 298947
# write_fst(AD_agg65, paste0(dir_data, "aggregated/AD_agg65_corrected.fst"))


##### Create first ADRD aggregated dataset, binning ages #####
# ADRD_agg <- dt[sex != 0 & race != 0 & ADRD_year >= cohort & ADRD_year >= year,
#              .(n_persons = .N, n_hosp = sum(ADRD_hosp), n_years = 1 + year - cohort),
#              by = .(zip, cohort, year, age_grp, sex, race, dual)]
# ADRD_agg[, .(.N, sum(n_hosp))] # 49531395 units 7480320 events
# write_fst(ADRD_agg, paste0(dir_proj, "data/aggregated/ADRD_agg_tv.fst")) # AD_agg_corrected.fst was pre-age-binning

# no unknown sex or race
# ADRD_agg <- ffs_entry_exit_adrd[sex != 0 & race != 0
#                                 , .(n_persons = .N, 
#                                     n_ADRDhosp = sum(ADRD_hosp), 
#                                     n_years = 1 + ADRD_year - ffs_entry_year),
#                                 by = .(ffs_entry_year, ADRD_zip, ADRD_age_binned, ADRD_year,
#                                        sexM, race_cat, any_dual)]
# ADRD_agg[, sum(n_ADRDhosp)] # 7540591
# write_fst(ADRD_agg, paste0(dir_proj, "data/aggregated/ADRD_agg_age_binned.fst")) # ADRD_agg_corrected.fst was pre-age-binning

# The following are equivalent to ADRD_agg (pre-age-binning)
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

# begin at age 65, no unknown sex or race (pre-age-binning)
# ADRD_agg65 <- ffs_entry_exit_adrd[entry_age == 65 & sex != 0 & race != 0
#                                 , .(n_persons = .N, 
#                                     n_ADRDhosp = sum(ADRD_hosp), 
#                                     n_years = 1 + ADRD_year - ffs_entry_year),
#                                 by = .(ffs_entry_year, ADRD_zip, ADRD_age, ADRD_year,
#                                        sexM, race_cat, any_dual)]
# ADRD_agg65[, sum(n_ADRDhosp)] # 968280
# write_fst(ADRD_agg65, paste0(dir_data, "aggregated/ADRD_agg65_corrected.fst"))
