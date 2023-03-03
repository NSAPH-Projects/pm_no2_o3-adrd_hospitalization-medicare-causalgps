# ------------------------------------------------------------------------------
#' Project: Causal ADRD
#' Code: determine FFS enrollment period for each individual
#' Inputs: denominator files, QD exposure data, meteorological data, hosp data
#' Outputs: qid entry and exit from ffs enrollment
#' Author: Daniel Mork
#' Last updated: May 12, 2022 (by Daniel Mork)
#' Converting into package format: March 02, 2022
#' Memory to run: 96 GB
# ------------------------------------------------------------------------------
#
# This process requires to have access to two external folders.
# ext_1 <- "/n/dominici_nsaph_l3/Lab/projects/pm_no2_o3-adrd_hosp-medicare-causalgps/data/"
# ext_2 <- "/n/dominici_nsaph_l3/Lab/projects/analytic/denom_by_year/"
#
#
# initiate the process ---------------------------------------------------------
tryCatch({
  # code that may raise an error or exception
  full_path <- whereami::thisfile_source()
  sp_name  <- sub("\\.R$", "", basename(full_path))
  print(paste("Project name (1): ", sp_name))
}, error = function(e) {
  # code to handle the error or exception
  file_name_env <- sys.frame(1)
  sp_name  <- sub("\\.R$", "", basename(file_name_env$fileName))
  print(paste("Project name (2): ", sp_name))
})

path_obj <- initialize_sub_project(sp_name = sp_name)

# Look at generated path -------------------------------------------------------
print(path_obj)

# Read denom files -------------------------------------------------------------
options(stringsAsFactors = FALSE)
data.table::setDTthreads(threads = 48)
fst::threads_fst(nr_of_threads = 48, reset_after_fork = TRUE)

cat("\n Reading denominator files...")
f <- list.files(file.path(path_obj$dir_data_private_ext_2, "denom_by_year"),
                 pattern = "\\.fst", full.names = TRUE)
myvars <- c("qid", "year", "zip", "hmo_mo", "age",
            "race", "sex", "dual", "dead", "fips")
dt <- data.table::rbindlist(lapply(f[2:18], fst::read_fst, columns = myvars,
                            as.data.table = TRUE))

data.table::setkey(dt, year, qid)
dt[,.N] # 651691916 total person years of data
# ! Check how many people removed
dt[year < 2000 | hmo_mo != 0, .N] # 162148450 person years removed

dt <- dt[year >= 2000 &
           hmo_mo == 0 & # retain FFS only (zero HMO months)
           race != 0 & sex != 0]
data.table::setkey(dt, qid, year)
dt[, hmo_mo := NULL]
dt <- unique(dt, by = c("qid", "year"))
cat(dt[,.N], "records\n") # 489543466

# Determine last year of continuous enrollment ---------------------------------

cat("\n Determine cohort and last year FFS")
# function to get last year of continuous enrollment in Medicare FFS
last_yr_cont <- function(year) {
  year[data.table::first(which(!(data.table::first(year):2017 %in% year))) - 1]
}
dt[, cohort := data.table::first(year), by = qid]
dt[, last_year_ffs := last_yr_cont(year), by = qid]
# remove records after year without some FFS coverage
dt[year > last_year_ffs, .N] # 12039782 removed
dt <- dt[year <= last_year_ffs]
dt[, .N] # 477503684
# corrected age (in case of inconsistencies)
dt[, age_corrected := data.table::first(age) + year - data.table::first(year),
   by = qid]
dt[, age := NULL]

# Merge RTI race code ----------------------------------------------------------
rti <- data.table::rbindlist(lapply(c(2009:2014, 2016), function(y) {
  d <- data.table::fread(file.path(path_obj$dir_data_private_ext_1,
                       "main_data_1_private",
                       "auxiliary_medicare_cols",
                       paste0("rti_race_", y, ".csv")))
  d[, year := y]
}), fill = TRUE)
rti[, rti_race_cd := as.integer(rti_race_cd)]
rti <- rti[!is.na(rti_race_cd) & rti_race_cd != 0]
dt <- merge(dt, rti, by = c("qid", "year"), all.x = TRUE)
dt[!is.na(rti_race_cd), race := rti_race_cd]
rm(rti)

# Finialize --------------------------------------------------------------------
finalize_subproject(path_obj = path_obj)



