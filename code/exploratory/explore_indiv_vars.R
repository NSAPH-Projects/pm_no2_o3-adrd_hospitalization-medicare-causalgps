#### Note: This script should be run after the data have been cleaned, processed, and aggregated data but NOT trimmed (i.e., files #1-5 in the code/aggregation folder have been run) ####


#### Set up computations ####

library(data.table)
library(fst)

# get directories and classifications of variables
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/git/code/"
source(paste0(dir_code, "constants.R"))

# optional: user sets the number of cores in this computing job
# more important: request enough memory (at least 100 GB) to open the large patient dataset
n_cores <- 48
setDTthreads(threads = n_cores)


### Get all patient data ###
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

# get data for each patient's FFS entry year only (i.e., ignoring later years of observation), for full cohort and ADRD cohort
dt_entry <- dt[year == cohort]
dt_ADRD_entry <- dt_entry[qid %in% ADRD_patients]

# collect relevant datasets into list for easy indexing
cohorts <- list(dt, dt_ADRD)
entry_year_cohorts <- list(dt_entry, dt_ADRD_entry)
names(cohorts) <- c("Full", "ADRD")
names(entry_year_cohorts) <- c("Full", "ADRD")


#### Set up txt file to store results ####
cat(paste("Variable",
          "Cohort",
          "Value",
          sep = ","),
    sep = "\n",
    file = paste0(dir_results, "exploratory/table1.txt"),
    append = TRUE)


#### Calculate and save number of individuals and person-years ####
for (cohort in c("Full", "ADRD")){
  cat(paste("Number of individuals",
            paste0(cohort, "Cohort"),
            uniqueN(cohorts[[cohort]][["qid"]]),
            sep = ","),
      sep = "\n",
      file = paste0(dir_results, "exploratory/table1.txt"),
      append = TRUE)
  cat(paste("Number of person-years",
            paste0(cohort, "Cohort"),
            nrow(cohorts[[cohort]]),
            sep = ","),
      sep = "\n",
      file = paste0(dir_results, "exploratory/table1.txt"),
      append = TRUE)
}


### Calculate and save individual-level variables, at each patient's year of FFS entry ###
for (cohort in c("Full", "ADRD")){
  
  # percent female (coded as 2 in the data) and male (coded as 1 in the data)
  cat(paste("Female",
            paste0(cohort, "Cohort"),
            paste0(round(mean(cohorts[[cohort]][["sex"]] == 2) * 100, 1), "%"),
            sep = ","),
      sep = "\n",
      file = paste0(dir_results, "exploratory/table1.txt"),
      append = TRUE)
  cat(paste("Male",
            paste0(cohort, "Cohort"),
            paste0(round(mean(cohorts[[cohort]][["sex"]] == 1) * 100, 1), "%"),
            sep = ","),
      sep = "\n",
      file = paste0(dir_results, "exploratory/table1.txt"),
      append = TRUE)
  
  # percent of each age group (5-year bins)
  cat(paste(levels(cohorts[[cohort]][["age_grp"]]),
            paste0(cohort, "Cohort"),
            paste0(round(prop.table(table(cohorts[[cohort]][["age_grp"]])) * 100, 1), "%"),
            sep = ","),
      sep = "\n",
      file = paste0(dir_results, "exploratory/table1.txt"),
      append = TRUE)
  
  # percent of each race (RTI-augmented race codes)
  cat(paste("Non-Hispanic White",
            paste0(cohort, "Cohort"),
            paste0(round(mean(cohorts[[cohort]][["race"]] == 1) * 100, 1), "%"),
            sep = ","),
      sep = "\n",
      file = paste0(dir_results, "exploratory/table1.txt"),
      append = TRUE)
  cat(paste("Black",
            paste0(cohort, "Cohort"),
            paste0(round(mean(cohorts[[cohort]][["race"]] == 2) * 100, 1), "%"),
            sep = ","),
      sep = "\n",
      file = paste0(dir_results, "exploratory/table1.txt"),
      append = TRUE)
  cat(paste("Hispanic",
            paste0(cohort, "Cohort"),
            paste0(round(mean(cohorts[[cohort]][["race"]] == 5) * 100, 1), "%"),
            sep = ","),
      sep = "\n",
      file = paste0(dir_results, "exploratory/table1.txt"),
      append = TRUE)
  cat(paste("Asian/Pacific Islander",
            paste0(cohort, "Cohort"),
            paste0(round(mean(cohorts[[cohort]][["race"]] == 4) * 100, 1), "%"),
            sep = ","),
      sep = "\n",
      file = paste0(dir_results, "exploratory/table1.txt"),
      append = TRUE)
  cat(paste("American Indian/Alaska Native",
            paste0(cohort, "Cohort"),
            paste0(round(mean(cohorts[[cohort]][["race"]] == 6) * 100, 1), "%"),
            sep = ","),
      sep = "\n",
      file = paste0(dir_results, "exploratory/table1.txt"),
      append = TRUE)
  cat(paste("Other",
            paste0(cohort, "Cohort"),
            paste0(round(mean(cohorts[[cohort]][["race"]] == 3) * 100, 1), "%"),
            sep = ","),
      sep = "\n",
      file = paste0(dir_results, "exploratory/table1.txt"),
      append = TRUE)
  
  ### to do: check and run the following
  
  # percent Medicaid eligible
  cat(paste("Not Medicaid-eligible",
            paste0(cohort, "Cohort"),
            paste0(round(mean(cohorts[[cohort]][["dual"]] == 0) * 100, 1), "%"),
            sep = ","),
      sep = "\n",
      file = paste0(dir_results, "exploratory/table1.txt"),
      append = TRUE)
  cat(paste("Medicaid eligible",
            paste0(cohort, "Cohort"),
            paste0(round(mean(cohorts[[cohort]][["dual"]] == 1) * 100, 1), "%"),
            sep = ","),
      sep = "\n",
      file = paste0(dir_results, "exploratory/table1.txt"),
      append = TRUE)
  
  ## to do: delete the following
  
  # prop.table(table(dt_entry$dual)) # Medicaid eligibility
}


### Calculate and save mean and SD of exposures and neighborhood variables, at each patient's year of FFS entry ###

# get exposures for each ZIP code and year
expos_dat <- read_fst(paste0(dir_data, "denom/year_zip_exposures.fst"),
                      as.data.table = TRUE)
expos_dat[, dat_year := year + 1][, year := NULL] # year for merging into dataset (year before patient enters FFS)
expos_dat[, zip := as.integer(zip)]
setkey(expos_dat, zip, dat_year)

# get confounders for each ZIP code and year
yr_zip_dat <- read_fst(paste0(dir_data, "denom/year_zip_confounders.fst"),
                       as.data.table = TRUE)
yr_zip_dat[, dat_year := year + 1][, year := NULL] # year for merging into dataset (year before patient enters FFS)
yr_zip_dat[, zip := as.integer(zip)]
setnames(yr_zip_dat,
         old = c("pct_blk", "pct_owner_occ"),
         new = c("prop_blk", "prop_owner_occ"))
setkey(yr_zip_dat, zip, dat_year)

# merge exposures with patient data, at each patient's FFS entry year
dt_entry <- merge(dt_entry, expos_dat,
                  by.x = c("zip", "year"), by.y = c("zip", "dat_year"),
                  all.x = TRUE)
dt_ADRD_entry <- merge(dt_ADRD_entry, expos_dat,
                       by.x = c("zip", "year"), by.y = c("zip", "dat_year"),
                       all.x = TRUE)

# merge confounders with patient data, at each patient's FFS entry year
dt_entry <- merge(dt_entry, yr_zip_dat,
                  by.x = c("zip", "year"), by.y = c("zip", "dat_year"),
                  all.x = TRUE)
dt_ADRD_entry <- merge(dt_ADRD_entry, yr_zip_dat,
                       by.x = c("zip", "year"), by.y = c("zip", "dat_year"),
                       all.x = TRUE)

# calculate and save mean and SD of exposures and neighborhood variables, at each patient's year of FFS entry
for (cohort in c("Full", "ADRD")){
  for (var in c(zip_expos_names, zip_quant_var_names)){
    cat(paste(var,
              paste0(cohort, "Cohort"),
              round(mean(cohorts[[cohort]][[var]], na.rm = T), 1), " (", round(sd(cohorts[[cohort]][[var]], na.rm = T), 1), ")",
              sep = ","),
        sep = "\n",
        file = paste0(dir_results, "exploratory/table1.txt"),
        append = TRUE)
  }
}


## to do: delete the following

# for (var in c(zip_expos_names, zip_quant_var_names)){
#   print(paste0(var, ": ", round(mean(dt_entry[[var]], na.rm = T), 1), " (", round(sd(dt_entry[[var]], na.rm = T), 1), ")"))
#   # cat("Mean of", var, "in full cohort at entry year:", mean(dt_entry[[var]], na.rm = T), "\n")
#   # cat("SD of", var, "in full cohort at entry year:", sd(dt_entry[[var]], na.rm = T), "\n")
# }
# 
# for (var in c(zip_expos_names, zip_quant_var_names)){
#   print(paste0(var, ": ", round(mean(dt_ADRD_entry[[var]], na.rm = T), 1), " (", round(sd(dt_ADRD_entry[[var]], na.rm = T), 1), ")"))
#   # cat("Mean of", var, "in ADRD cohort at entry year:", mean(dt_ADRD_entry[[var]], na.rm = T), "\n")
#   # cat("SD of", var, "in ADRD cohort at entry year:", sd(dt_ADRD_entry[[var]], na.rm = T), "\n")
# }

