### Note: This script should be run after the data have been cleaned, processed, and aggregated data but NOT trimmed (i.e., files #1-5 in the code/aggregation folder have been run) ###

library(data.table)
library(fst)

# get directories and classifications of variables
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/git/code/"
source(paste0(dir_code, "constants.R"))

# user should set number of cores in this computing job
n_cores <- 48
setDTthreads(threads = n_cores)


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

# get subset of patients who experienced ADRD event
dt_ADRD <- dt[year == ADRD_year][ADRD_hosp == 1]
ADRD_patients <- dt_ADRD$qid
dt_ADRD <- dt[qid %in% ADRD_patients]

# get patient-level variables at each patient's entry year, for full cohort and ADRD cohort
dt_entry <- dt[year == cohort]
dt_ADRD_entry <- dt_entry[qid %in% ADRD_patients]

# collect relevant datasets into list for easy indexing
cohorts <- list(dt, dt_ADRD)
entry_year_cohorts <- list(dt_entry, dt_ADRD_entry)
names(cohorts) <- c("Full", "ADRD")
names(entry_year_cohorts) <- c("Full", "ADRD")

# set up txt file to store results
cat(paste("Variable", "Cohort", "Value", sep = ","),
    sep = "\n",
    file = paste0(dir_results, "exploratory/table1.txt"),
    append = TRUE)

# save number of individuals and person-years
for (cohort in c("Full", "ADRD")){
  cat(paste("Number of individuals", paste0(cohort, "Cohort"), uniqueN(cohorts[[cohort]][["qid"]]), sep = ","),
      sep = "\n",
      file = paste0(dir_results, "exploratory/table1.txt"),
      append = TRUE)
  cat(paste("Number of person-years", paste0(cohort, "Cohort"), nrow(cohorts[[cohort]]), sep = ","),
      sep = "\n",
      file = paste0(dir_results, "exploratory/table1.txt"),
      append = TRUE)
}

# save individual-level variables at each patient's year of entry
for (cohort in c("Full", "ADRD")){
  
  # proportion female (coded as 2) and male (coded as 1)
  cat(paste("Female", paste0(cohort, "Cohort"), paste0(round(mean(cohorts[[cohort]][["sex"]] == 2) * 100, 1), "%"), sep = ","),
      sep = "\n",
      file = paste0(dir_results, "exploratory/table1.txt"),
      append = TRUE)
  cat(paste("Male", paste0(cohort, "Cohort"), paste0(round(mean(cohorts[[cohort]][["sex"]] == 1) * 100, 1), "%"), sep = ","),
      sep = "\n",
      file = paste0(dir_results, "exploratory/table1.txt"),
      append = TRUE)
  
  # proportion of each age group (5-year bins)
  cat(paste(levels(cohorts[[cohort]][["age_grp"]]), paste0(cohort, "Cohort"), paste0(round(prop.table(table(cohorts[[cohort]][["age_grp"]])) * 100, 1), "%"), sep = ","),
      sep = "\n",
      file = paste0(dir_results, "exploratory/table1.txt"),
      append = TRUE)
  
  # RTI-augmented race codes
  cat(paste("Non-Hispanic White", paste0(cohort, "Cohort"), paste0(round(mean(cohorts[[cohort]][["race"]] == 1) * 100, 1), "%"), sep = ","),
      sep = "\n",
      file = paste0(dir_results, "exploratory/table1.txt"),
      append = TRUE)
  cat(paste("Black", paste0(cohort, "Cohort"), paste0(round(mean(cohorts[[cohort]][["race"]] == 2) * 100, 1), "%"), sep = ","),
      sep = "\n",
      file = paste0(dir_results, "exploratory/table1.txt"),
      append = TRUE)
  cat(paste("Hispanic", paste0(cohort, "Cohort"), paste0(round(mean(cohorts[[cohort]][["race"]] == 5) * 100, 1), "%"), sep = ","),
      sep = "\n",
      file = paste0(dir_results, "exploratory/table1.txt"),
      append = TRUE)
  cat(paste("Asian/Pacific Islander", paste0(cohort, "Cohort"), paste0(round(mean(cohorts[[cohort]][["race"]] == 4) * 100, 1), "%"), sep = ","),
      sep = "\n",
      file = paste0(dir_results, "exploratory/table1.txt"),
      append = TRUE)
  cat(paste("American Indian/Alaska Native", paste0(cohort, "Cohort"), paste0(round(mean(cohorts[[cohort]][["race"]] == 6) * 100, 1), "%"), sep = ","),
      sep = "\n",
      file = paste0(dir_results, "exploratory/table1.txt"),
      append = TRUE)
  cat(paste("Other", paste0(cohort, "Cohort"), paste0(round(mean(cohorts[[cohort]][["race"]] == 3) * 100, 1), "%"), sep = ","),
      sep = "\n",
      file = paste0(dir_results, "exploratory/table1.txt"),
      append = TRUE)
}

## to do: everything below this


prop.table(table(dt_entry$dual)) # Medicaid eligibility



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

for (var in c(zip_expos_names, zip_quant_var_names)){
  print(paste0(var, ": ", round(mean(dt_ADRD_entry[[var]], na.rm = T), 1), " (", round(sd(dt_ADRD_entry[[var]], na.rm = T), 1), ")"))
  # cat("Mean of", var, "in ADRD cohort at entry year:", mean(dt_ADRD_entry[[var]], na.rm = T), "\n")
  # cat("SD of", var, "in ADRD cohort at entry year:", sd(dt_ADRD_entry[[var]], na.rm = T), "\n")
}

### Not used in manuscript

round(prop.table(table(dt_entry$region)), 3)
round(prop.table(table(dt_ADRD_entry$region)), 3)

# get distribution of patients' years of followup
followup_years <- dt[, .N, by = qid]
summary(followup_years$N) # min is 1, max is 16, median is 7, mean is 8.296

# for patients who experienced ADRD event, get distribution of patients' years of followup
ADRD_patient_followup_years <- dt_ADRD[, .N, by = qid]
summary(ADRD_patient_followup_years$N) # min is 2, max is 16, median is 7, mean is 7.933