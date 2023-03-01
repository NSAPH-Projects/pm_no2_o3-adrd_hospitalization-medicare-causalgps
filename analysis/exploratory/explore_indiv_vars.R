library(data.table)
library(fst)

# directories for data, code, and results
dir_data <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/data/"
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/code/"
dir_results <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/results/"

# user should set number of cores in this computing job
n_cores <- 16
setDTthreads(threads = n_cores)

# get classifications of variables
source(paste0(dir_code, "analysis/helper_functions.R"))

# read in raw (time-varying) data
dt <- read_fst(paste0(dir_data, "denom/complete_ADRD_denom.fst"),
               as.data.table = TRUE)
setkey(dt, zip, cohort, year, age_grp, sex, race, dual)

# get (time-varying) data for our analyses
dt <- dt[ADRD_year >= cohort + 2 & # 2 year 'clean' period
           ADRD_year >= year & # censor observations after ADRD event
           year >= 2001] # begin with year 2001 to pair with 2000 exposures
dt[cohort == 2000, cohort := 2001] # start 2000 cohort in 2001

# calculate number of individuals, person-years, events
cat("Number of individuals:", uniqueN(dt$qid)) # 50,053,399
cat("Number of person-years:", nrow(dt)) # 415,260,761
cat("Number of events:", sum(dt$ADRD_hosp)) # 5,935,558
cat("Overall prevalence across person-years:", sum(dt$ADRD_hosp) / nrow(dt)) # 0.01429357

# calculate median years of followup
# to do. need to figure out which of these data.frames to use
n_years_followup1 <- dt[year == ADRD_year] # 49,589,968 rows. unique qid's. year of event or censoring
n_years_followup2 <- dt[year == max(year), .SD, by = qid] # 22,490,481 rows. unique qid's. should be equivalent to #1?
# n_years_followup3 <- dt[year == last_year_ffs] # 45,995,856 rows. unique qid's. year of exiting ffs, which may be after event or censoring

# to do. need to figure out which of these data.frames to use
n_years_followup1[, n_years := year - cohort + 1]
summary(n_years_followup1$n_years)

# alternatively (equivalent?)
dt[, .N, by = qid] # will give the number of followup years by individual
dt[, .N, by = qid][, median(N)] # should give the median follow up years


# calculate summary statistics of patient-level variables at each patient's entry year
# note: this is for the untrimmed exposure
dt_entry <- dt[year == cohort] # each patient's entry year

table(dt_entry$sex) # 1 = male, 2 = female
table(dt_entry$age_grp) # in 5-year bins
table(dt_entry$race) # RTI-augmented race codes: https://resdac.org/cms-data/variables/research-triangle-institute-rti-race-code
table(dt_entry$dual) # Medicaid eligibility

prop.table(table(dt_entry$sex)) # 1 = male, 2 = female
prop.table(table(dt_entry$age_grp)) # in 5-year bins
prop.table(table(dt_entry$race)) # RTI-augmented race codes: https://resdac.org/cms-data/variables/research-triangle-institute-rti-race-code
prop.table(table(dt_entry$dual)) # Medicaid eligibility