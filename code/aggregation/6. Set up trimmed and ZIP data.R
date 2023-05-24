rm(list = ls())
gc()

##### 1. Setup #####
library(data.table)
library(fst)

# directories for data, code, and results
dir_data <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/data/"
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/code/"
dir_results <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/results/"

# read in full data
ADRD_agg_lagged <- read_fst(paste0(dir_data, "analysis/ADRD_complete_tv.fst"),
                            as.data.table = TRUE)
n_rows_untrimmed <- nrow(ADRD_agg_lagged) # 34,763,397 rows
setnames(ADRD_agg_lagged,
         old = c("pct_blk", "pct_owner_occ"),
         new = c("prop_blk", "prop_owner_occ"))
ADRD_agg_lagged[, sex := ifelse(sex == 2, 0, # in raw data, female is coded as 2; recode as 0
                               ifelse(sex == 1, 1, NA))] # in raw data, male is coded as 1; keep this; there should be no NAs
ADRD_agg_lagged[, `:=`(zip = as.factor(zip),
                       year = as.factor(year),
                       cohort = as.factor(cohort),
                       age_grp = as.factor(age_grp),
                       sex = as.factor(sex),
                       race = as.factor(race),
                       dual = as.factor(dual))]

source(paste0(dir_code, "analysis/helper_functions.R"))


##### 2. Get data for exposure, outcome, and covariates of interest #####

exposure_name <- "pm25"

ADRD_agg_lagged_subset <- subset(ADRD_agg_lagged, select = c(exposure_name,
                                                             outcome_name,
                                                             zip_var_names,
                                                             "zip",
                                                             indiv_var_names,
                                                             offset_var_names))
for (var in c("zip",
              zip_unordered_cat_var_names,
              indiv_unordered_cat_var_names)){
  ADRD_agg_lagged_subset[[var]] <- as.factor(ADRD_agg_lagged_subset[[var]])
}

# rename exposure column to "w" and outcome column to "Y"
colnames(ADRD_agg_lagged_subset)[colnames(ADRD_agg_lagged_subset) == exposure_name] <- "w"
colnames(ADRD_agg_lagged_subset)[colnames(ADRD_agg_lagged_subset) == outcome_name] <- "Y"

# for (var in zip_ordered_cat_var_names){
#   ADRD_agg_lagged_subset[[var]] <- factor(ADRD_agg_lagged_subset[[var]], ordered = TRUE)
# }


##### 3. Isolate just the exposure and covariates (1 observation per ZIP, year) #####

zip_year_data <- subset(ADRD_agg_lagged_subset, select = c("zip",
                                                           "year",
                                                           "w",
                                                           zip_var_names))
zip_year_data <- unique(zip_year_data, by = c("zip", "year"))
n_zip_year_rows_untrimmed <- nrow(zip_year_data) # 486,793 rows

##### 4. Trim exposure #####

trim_level <- c(0.05, 0.95)
save_fst <- T

trim <- quantile(zip_year_data$w, trim_level)
zip_year_rows_within_range <- zip_year_data$w >= trim[1] & zip_year_data$w <= trim[2]
zip_year_trimmed <- zip_year_data[zip_year_rows_within_range, ]
n_zip_year_rows_trimmed <- nrow(zip_year_trimmed) # 438,113 rows

if (save_fst){
  write_fst(zip_year_trimmed, paste0(dir_data,
                                     "analysis/",
                                     exposure_name, "/",
                                     "zip_year_data_trimmed_", trim_level[1], "_", trim_level[2],
                                     ".fst")) 
}

ADRD_agg_rows_within_range <- ADRD_agg_lagged_subset$w >= trim[1] & ADRD_agg_lagged_subset$w <= trim[2]
ADRD_agg_lagged_trimmed <- ADRD_agg_lagged_subset[ADRD_agg_rows_within_range, ]
n_rows_trimmed <- nrow(ADRD_agg_lagged_trimmed) # 31,989,354 for PM2.5; 30,813507 for NO2

if (save_fst){
  write_fst(ADRD_agg_lagged_trimmed, paste0(dir_data,
                                            "analysis/",
                                            exposure_name, "/",
                                            "zip_year_data_with_strata_trimmed_", trim_level[1], "_", trim_level[2],
                                            ".fst"))
}