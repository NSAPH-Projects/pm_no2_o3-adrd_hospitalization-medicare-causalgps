###############################################################################
# Project: Causal effect of air pollution on first ADRD hospitalization       #
# Code: extract follow up period using denominator files                      #
# Input: "First_hosp_yr.fst"                                                  #
# Input: denominator files                                                    #                           
# Output: "ADRD_aggregated_entry_to_followup.fst"                             #
# Output: "zip_year_census_data.fst"                                          #
# Author: Daniel Mork                                                         #
# Date: 2021-12-21                                                            #
###############################################################################
rm(list = ls())
gc()
cores = 4
############################# 0. Setup ########################################
library(data.table)
library(fst)
library(dplyr)
library(NSAPHutils)
options(stringsAsFactors = FALSE)

setDTthreads(threads = cores)
dir_data <- "/nfs/home/D/dam9096/shared_space/ci3_analysis/dmork/Data/Causal_ADRD/"
dir_denominator <- "/nfs/home/D/dam9096/shared_space/ci3_health_data/medicare/mortality/1999_2016/wu/cache_data/merged_by_year_v2/"

f <- list.files(dir_denominator, pattern = "\\.fst", full.names = TRUE)
# example <- read_fst(f[1], from = 1, to = 10000, as.data.table = TRUE)
# names(example)
# example
# > names(example)
# [1] "zip"                          "year"                         "qid"                          "dodflag"                     
# [5] "bene_dod"                     "sex"                          "race"                         "age"                         
# [9] "hmo_mo"                       "hmoind"                       "statecode"                    "latitude"                    
# [13] "longitude"                    "dual"                         "death"                        "dead"                        
# [17] "entry_age"                    "entry_year"                   "entry_age_break"              "followup_year"               
# [21] "followup_year_plus_one"       "pm25_ensemble"                "pm25_no_interp"               "pm25_nn"                     
# [25] "ozone"                        "ozone_no_interp"              "zcta"                         "poverty"                     
# [29] "popdensity"                   "medianhousevalue"             "pct_blk"                      "medhouseholdincome"          
# [33] "pct_owner_occ"                "hispanic"                     "education"                    "population"                  
# [37] "zcta_no_interp"               "poverty_no_interp"            "popdensity_no_interp"         "medianhousevalue_no_interp"  
# [41] "pct_blk_no_interp"            "medhouseholdincome_no_interp" "pct_owner_occ_no_interp"      "hispanic_no_interp"          
# [45] "education_no_interp"          "population_no_interp"         "smoke_rate"                   "mean_bmi"                    
# [49] "smoke_rate_no_interp"         "mean_bmi_no_interp"           "amb_visit_pct"                "a1c_exm_pct"                 
# [53] "amb_visit_pct_no_interp"      "a1c_exm_pct_no_interp"        "tmmx"                         "rmax"                        
# [57] "pr"                           "cluster_cat"                  "fips_no_interp"               "fips"                        
# [61] "summer_tmmx"                  "summer_rmax"                  "winter_tmmx"                  "winter_rmax" 
myvars <- c("qid", "year", "zip", "sex", "race", "age", "dual")
qid_dat <- rbindlist(lapply(f, function(x) {
  read_fst(x, columns = myvars, as.data.table = TRUE)
}, mc.cores = cores))
write_fst(qid_dat, paste0(dir_data, "Denominator/qid_complete.fst"))







myvars <- c("qid", "year", "zip", "sex", "race", "age", "dual", "statecode", 
            "dead", "mean_bmi", "smoke_rate", "hispanic", "pct_blk", 
            "medhouseholdincome", "medianhousevalue", "poverty", "education", 
            "popdensity", "pct_owner_occ", "summer_tmmx", "summer_rmax",
            "winter_tmmx", "winter_rmax")
state48 <- c("AL", "AZ", "AR", "CA", "CO", "CT", "DE", "DC", "FL", "GA", #includes Wash DC
             "ID", "IL", "IN", "IA", "KS", "KY", "LA", "ME", "MD", "MA", "MI",
             "MN", "MS", "MO", "MT", "NE", "NV", "NH", "NJ", "NM", "NY", "NC",
             "ND", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN", "TX", "UT",
             "VT", "VA", "WA", "WV", "WI", "WY")

ADRDhosp <- read_fst(paste0(dir_data, "First_hosp_yr.fst"), as.data.table = TRUE)
dim(ADRDhosp)[1]
# [1] 7,249,590

########################### 1. Load denominator ################################
denomDat <- data.table(qid = character(), 
                       entryYear = integer(), followupYear = integer(), 
                       age = integer(), sex = integer(), race = character(),
                       dual = integer(), 
                       followupAge = integer(), followupDual = integer(),
                       entryZip = character(), followupZip = character(),
                       ADRDevent = logical(), ADRDyear = integer(),
                       dead = logical(), deathYear = integer())
zipYrDenom <- list()
# read in denominator files
for (i in 2:length(f)) { # begin with year 2000
  cat("Reading denominator file", f[i])
  dt <- read_fst(f[i], columns = myvars, as.data.table = TRUE)
  setnames(dt, "year", "entryYear")
  dt[, followupYear := entryYear]
  cat(" ... checkpoint: 1")
  
  # Remove impossible zips (< 00501), limit to contiguous 48 states, race != 0,
  # current year only, no people with previous hospitalization or death and duplicates
  # Recode race to wht, blk, his, oth, create 5-digit zipcode
  dt <- unique(dt[zip >= 501 & race != 0 & statecode %in% state48 & 
                    entryYear == 1998 + i &
                    !(qid %in% ADRDhosp[ADRDhosp$firstADRDyr < 1998 + i, QID]) &
                    !(qid %in% denomDat[deathYear < 1998 + i, qid])][, 
                  curzip := int_to_zip_str(zip)][,
                  race2 := ifelse(race == 1, "wht",
                                  ifelse(race == 2, "blk",
                                         ifelse(race == 5, "his", "oth")))][, 
                  race := NULL])
  dt[sex == 2, sex := 0] # recode 0 = female, 1 = male
  dt$dual <- as.integer(dt$dual) # recode dual to integer
  setnames(dt, c("zip", "race2"), c(paste0("zip", 1998 + i), "race"))
  cat(" 2")
  
  # Add ADRD event and year, death year to new data
  dt[, ADRDevent := FALSE]
  dt[qid %in% ADRDhosp[firstADRDyr == 1998 + i, QID], ADRDevent := TRUE]
  dt[ADRDevent == TRUE, ADRDyear := 1998 + i]
  dt[dead == TRUE, deathYear := 1998 + i]
  cat(" 3")

  # Update previous entries
  # update followup year if individual still in medicare
  denomDat[dead == FALSE & qid %in% dt$qid, followupYear := 1998 + i]
  # update ADRDevent if recorded hopsitalization
  denomDat[qid %in% ADRDhosp[firstADRDyr == 1998 + i, QID], ADRDevent := TRUE]
  # update ADRDyear if added event
  denomDat[ADRDevent == TRUE & is.na(ADRDyear), ADRDyear := 1998 + i]
  # update death of recorded in current year's data
  denomDat[qid %in% dt[dead == TRUE, qid], dead := TRUE]
  # update year of death
  denomDat[dead == TRUE & is.na(deathYear), deathYear := 1998 + i]
  cat(" 4")

  # Append new QIDs onto denomDat
  denomDat <- rbindlist(list(denomDat,
                             dt[!(qid %in% denomDat[, qid]), # new QIDs only
                                c("qid", "entryYear", "followupYear", "age", 
                                  "sex", "race", "dual", 
                                  "ADRDevent", "ADRDyear", 
                                  "dead", "deathYear"),
                                with = FALSE]),
                        use.names = TRUE, fill = TRUE)
  cat(" 5")

  # Merge current zips into denominator data
  setnames(dt, c("age", "dual"), c("age2", "dual2"))
  denomDat <- merge(denomDat,
                    dt[, c("qid", "curzip", "age2", "dual2", paste0("zip", 1998 + i)), 
                       with = FALSE],
                    by = "qid", all.x = TRUE)
  # set zipcode at entry into dataset
  denomDat[is.na(entryZip), entryZip := curzip]
  # set zipcode at end of followup (ADRD event or death)
  denomDat[is.na(followupZip) & !is.na(curzip), followupZip := curzip]
  # update followup age and dual
  denomDat[is.na(followupAge), followupAge := age2]
  denomDat[is.na(followupDual), followupDual := dual2]
  denomDat[, age2 := NULL][, dual2 := NULL][, curzip := NULL]
  cat(" 6")
  
  # Save aggregated zip/year data
  setnames(dt, c("age2", "dual2"), c("age", "dual"))
  zipYrDenom[[i]] <- dt[,.(.N,
                           "avg_age" = mean(age),
                           "med_age" = median(age),
                           "bmi" = mean(mean_bmi, na.rm = T),
                           "smoke" = mean(smoke_rate, na.rm = T),
                           "hispanic" = mean(hispanic, na.rm = T),
                           "pct_blk" = mean(pct_blk, na.rm = T),
                           "income" = mean(medhouseholdincome, na.rm = T),
                           "value" = mean(medianhousevalue, na.rm = T),
                           "poverty" = mean(poverty, na.rm = T),
                           "education" = mean(education, na.rm = T),
                           "popdensity" = mean(popdensity, na.rm = T),
                           "ownerocc" = mean(pct_owner_occ, na.rm = T),
                           "summer_tmmx" = mean(summer_tmmx, na.rm = T),
                           "summer_rmax" = mean(summer_rmax, na.rm = T),
                           "winter_tmmx" = mean(winter_tmmx, na.rm = T),
                           "winter_rmax" = mean(winter_rmax, na.rm = T),
                           "death_prob" = mean(dead, na.rm = T)),
                        by = c("entryYear", "curzip", "sex", "race", "dual")]
  cat(" done!\n")
  rm(dt)
  gc()
}

# set follow up zip to last year of recorded data
denomDat[is.na(followupZip), 
         followupZip := 
          int_to_zip_str(ifelse(!is.na(zip2016), zip2016,
                          ifelse(!is.na(zip2015), zip2015,
                           ifelse(!is.na(zip2014), zip2014,
                            ifelse(!is.na(zip2013), zip2013,
                             ifelse(!is.na(zip2012), zip2012,
                              ifelse(!is.na(zip2011), zip2011,
                               ifelse(!is.na(zip2010), zip2010,
                                ifelse(!is.na(zip2009), zip2009,
                                 ifelse(!is.na(zip2008), zip2008,
                                  ifelse(!is.na(zip2007), zip2007,
                                   ifelse(!is.na(zip2006), zip2006,
                                    ifelse(!is.na(zip2005), zip2005,
                                     ifelse(!is.na(zip2004), zip2004,
                                      ifelse(!is.na(zip2003), zip2003,
                                       ifelse(!is.na(zip2002), zip2002,
                                        ifelse(!is.na(zip2001), zip2001, zip2000)
                                        ))))))))))))))))]
                                     



############################# 2. Summarize denominator #########################
# Aggregated year/zip data
zipYrDenom <- rbindlist(zipYrDenom)
sum(zipYrDenom$N)
# [1] 615,216,264
setnames(zipYrDenom, "entryYear", "year")
zipYrDenom[, zip := int_to_zip_str(as.integer(curzip))][, curzip := NULL]
write_fst(zipYrDenom, paste0(dir_data, "zip_year_census_data.fst"))

# Complete data from denominator and hospitalizations (including possible missing)
denomDat[, .(.N, sum(followupYear - entryYear + 1))]
# 617942011 person-years of data
# 73759420 unique individuals
denomDat[ADRDevent == T, .(.N, mean(followupYear - entryYear) + 1,
                           mean(age + followupYear - entryYear))]
# 7219250 first-time ADRD events (9.7% of individuals)
# 7.166948 average years of follow-up to ADRD event
# 82.62138 mean age of ADRD event (year prior to event)
denomDat[dead == T, .(.N, mean(age + followupYear - entryYear))]
# 25868796 death events (1/3 of individuals)
# 80.86726 mean age of death (*note this is before mean age of ADRD event)
denomDat[dead == T & ADRDevent == T, .(.N, mean(age + followupYear - entryYear))]
# 11902223 individuals with ADRD event in year of death
# 84.54565 mean age of ADRD event and death in same year (year prior to event)
denomDat[dead == T & ADRDevent == F, .(.N, mean(age + followupYear - entryYear))]
# 23266292 individuals died without any recorded ADRD events
# 80.47863 mean age of death (year prior to event)
write_fst(denomDat, paste0(dir_data, "ADRD_complete_entry_to_followup.fst"))



################ 3. Aggregated denominator data for ADRD events ################
ADRD_aggregated <- denomDat[, .(nADRDEvents = sum(ADRDevent), nDeaths = sum(dead),
                                .N, nYears = followupYear - entryYear + 1),
                            by = c("entryYear", "followupYear", "followupZip",
                                   "age", "sex", "race", "dual")]
ADRD_aggregated[, sum(nYears * N)]
# [1] 615,439,100 # larger because some individuals not continuously enrolled
ADRD_aggregated[, followupAge := age + (followupYear - entryYear)]
write_fst(ADRD_aggregated, paste0(dir_data, "ADRD_aggregated_entry_to_followup.fst"))
