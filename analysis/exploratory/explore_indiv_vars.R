library(data.table)
library(fst)

# directories for data, code, and results
dir_data <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/data/"
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/code/"
dir_results <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/results/"

# user should set number of cores in this computing job
n_cores <- 48
setDTthreads(threads = n_cores)

# get classifications of variables
source(paste0(dir_code, "analysis/helper_functions.R"))

# get (time-varying) patient data
dt <- read_fst(paste0(dir_data, "denom/complete_ADRD_denom.fst"),
               as.data.table = TRUE)
setkey(dt, zip, cohort, year, age_grp, sex, race, dual)
dt <- dt[ADRD_year >= cohort + 2 & # 2 year 'clean' period
           ADRD_year >= year & # censor observations after ADRD event
           year >= 2001] # begin with year 2001 to pair with 2000 exposures
dt[cohort == 2000, cohort := 2001] # start 2000 cohort in 2001
dt <- subset(dt, select = c("qid",
                            "cohort",
                            "ADRD_year",
                            "last_year_ffs",
                            "ADRD_hosp",
                            "zip",
                            indiv_var_names)) # includes "year"

## Note: No trimming (by exposure or GPS) has been performed at this point

# calculate number of individuals, person-years, events
cat("Number of individuals:", uniqueN(dt$qid)) # 50,053,399
cat("Number of person-years:", nrow(dt)) # 415,260,761
cat("Number of events:", sum(dt$ADRD_hosp)) # 5,935,558
cat("Overall prevalence across person-years:", sum(dt$ADRD_hosp) / nrow(dt)) # 0.01429357

# get distribution of patients' years of followup
followup_years <- dt[, .N, by = qid]
summary(followup_years$N) # min is 1, max is 16, median is 7, mean is 8.296

# calculate summary statistics of patient-level variables at each patient's entry year
dt_entry <- dt[year == cohort]
prop.table(table(dt_entry$sex))
prop.table(table(dt_entry$age_grp)) # in 5-year bins
prop.table(table(dt_entry$race)) # RTI-augmented race codes
prop.table(table(dt_entry$dual)) # Medicaid eligibility

# get subset of patients who experienced ADRD event
dt_ADRD <- dt[year == ADRD_year][ADRD_hosp == 1]
ADRD_patients <- dt_ADRD$qid
dt_ADRD <- dt[qid %in% ADRD_patients]

# for patients who experienced ADRD event, calculate number of individuals, person-years, events
cat("Number of person-years among patients with ADRD event:", nrow(dt_ADRD)) # 47,084,973
cat("Overall prevalence across person-years:", sum(dt$ADRD_hosp) / nrow(dt_ADRD)) # 0.01429357

# for patients who experienced ADRD event, get distribution of patients' years of followup
ADRD_patient_followup_years <- dt_ADRD[, .N, by = qid]
summary(ADRD_patient_followup_years$N) # min is 2, max is 16, median is 7, mean is 7.933

# for patients who experienced ADRD event, calculate summary statistics of patient-level variables at each patient's entry year
dt_ADRD_entry <- dt_entry[qid %in% ADRD_patients]
prop.table(table(dt_ADRD_entry$sex))
prop.table(table(dt_ADRD_entry$age_grp)) # in 5-year bins
prop.table(table(dt_ADRD_entry$race)) # RTI-augmented race codes
prop.table(table(dt_ADRD_entry$dual)) # Medicaid eligibility

#### Merge exposures and confounders, at entry year ####

yr_zip_dat <- read_fst(paste0(dir_data, "denom/year_zip_confounders.fst"),
                       as.data.table = TRUE)
yr_zip_dat[, dat_year := year + 1][, year := NULL] # year for merging into dataset (year before end)
yr_zip_dat[, zip := as.integer(zip)]
setnames(yr_zip_dat,
         old = c("pct_blk", "pct_owner_occ"),
         new = c("prop_blk", "prop_owner_occ"))
setkey(yr_zip_dat, zip, dat_year)

expos_dat <- read_fst(paste0(dir_data, "denom/year_zip_exposures.fst"),
                      as.data.table = TRUE)
expos_dat[, dat_year := year + 1][, year := NULL] # year for merging into dataset (year before end)
expos_dat[, zip := as.integer(zip)]
setkey(expos_dat, zip, dat_year)

dt_entry <- merge(dt_entry, yr_zip_dat,
                  by.x = c("zip", "year"), by.y = c("zip", "dat_year"),
                  all.x = TRUE)
dt_entry <- merge(dt_entry, expos_dat,
                  by.x = c("zip", "year"), by.y = c("zip", "dat_year"),
                  all.x = TRUE)

dt_ADRD_entry <- merge(dt_ADRD_entry, yr_zip_dat,
                       by.x = c("zip", "year"), by.y = c("zip", "dat_year"),
                       all.x = TRUE)
dt_ADRD_entry <- merge(dt_ADRD_entry, expos_dat,
                       by.x = c("zip", "year"), by.y = c("zip", "dat_year"),
                       all.x = TRUE)

for (var in c(zip_expos_names, zip_quant_var_names)){
  print(paste0(var, ": ", round(mean(dt_entry[[var]], na.rm = T), 1), " (", round(sd(dt_entry[[var]], na.rm = T), 1), ")"))
  # cat("Mean of", var, "in full cohort at entry year:", mean(dt_entry[[var]], na.rm = T), "\n")
  # cat("SD of", var, "in full cohort at entry year:", sd(dt_entry[[var]], na.rm = T), "\n")
}

round(prop.table(table(dt_entry$region)), 3)

for (var in c(zip_expos_names, zip_quant_var_names)){
  print(paste0(var, ": ", round(mean(dt_ADRD_entry[[var]], na.rm = T), 1), " (", round(sd(dt_ADRD_entry[[var]], na.rm = T), 1), ")"))
  # cat("Mean of", var, "in ADRD cohort at entry year:", mean(dt_ADRD_entry[[var]], na.rm = T), "\n")
  # cat("SD of", var, "in ADRD cohort at entry year:", sd(dt_ADRD_entry[[var]], na.rm = T), "\n")
}

round(prop.table(table(dt_ADRD_entry$region)), 3)
