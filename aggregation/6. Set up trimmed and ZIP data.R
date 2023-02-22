rm(list = ls())
gc()

##### 0. Setup #####
# devtools::install_github("fasrc/CausalGPS", ref="develop")
library(data.table)
library(fst)

# directories for data, code, and results
dir_data <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/data/"
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/code/"
dir_results <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/results/"

# read in full data
ADRD_agg_lagged <- read_fst(paste0(dir_data, "analysis/ADRD_complete_tv.fst"), as.data.table = TRUE)
n_rows_untrimmed <- nrow(ADRD_agg_lagged) # 34,763,397 rows
setnames(ADRD_agg_lagged,
         old = c("pct_blk", "pct_owner_occ"),
         new = c("prop_blk", "prop_owner_occ"))
ADRD_agg_lagged[, `:=`(zip = as.factor(zip),
                       year = as.factor(year),
                       cohort = as.factor(cohort),
                       age_grp = as.factor(age_grp),
                       sex = as.factor(sex),
                       race = as.factor(race),
                       dual = as.factor(dual))]

source(paste0(dir_code, "analysis/helper_functions.R"))


##### Get data for exposure, outcome, and covariates of interest #####

exposure_name <- "pm25"
other_expos_names <- zip_expos_names[zip_expos_names != exposure_name]
outcome_name <- "n_hosp"
ADRD_agg_lagged_subset <- subset(ADRD_agg_lagged,
                                             select = c(exposure_name, outcome_name, other_expos_names, zip_var_names, "zip", indiv_var_names, offset_var_names))
for (var in c("zip", zip_unordered_cat_var_names, indiv_unordered_cat_var_names)){
  ADRD_agg_lagged_subset[[var]] <- as.factor(ADRD_agg_lagged_subset[[var]])
}
colnames(ADRD_agg_lagged_subset)[colnames(ADRD_agg_lagged_subset) == exposure_name] <- "w"
# for (var in zip_ordered_cat_var_names){
#   ADRD_agg_lagged_subset[[var]] <- factor(ADRD_agg_lagged_subset[[var]], ordered = TRUE)
# }


##### Isolate just the exposure and covariates (1 observation per ZIP, year) and trim exposure #####

zip_year_data <- subset(ADRD_agg_lagged_subset,
                                    select = c("zip", "year", "w", other_expos_names, zip_var_names))
zip_year_data <- unique(zip_year_data, by = c("zip", "year"))
n_zip_year_rows_untrimmed <- nrow(zip_year_data) # 486,793 rows

trim_1_99 <- quantile(zip_year_data$w, c(0.01, 0.99))
zip_year_rows_within_range <- zip_year_data$w >= trim_1_99[1] & zip_year_data$w <= trim_1_99[2]
zip_year_trimmed_1_99 <- zip_year_data[zip_year_rows_within_range, ]
n_zip_year_rows <- nrow(zip_year_trimmed_1_99)
# write_fst(zip_year_trimmed_1_99, paste0(dir_data, "analysis/pm_zip_year_data_trimmed_1_99.fst"))

ADRD_agg_rows_within_range <- ADRD_agg_lagged_subset$w >= trim_1_99[1] & ADRD_agg_lagged_subset$w <= trim_1_99[2]
ADRD_agg_lagged_trimmed_1_99 <- ADRD_agg_lagged_subset[ADRD_agg_rows_within_range, ]
n_rows <- nrow(ADRD_agg_lagged_trimmed_1_99) # 34,141,155 for PM1.5; 34,090,022 for NO2; 33,411,373 for ozone
# write_fst(ADRD_agg_lagged_trimmed_1_99, paste0(dir_data, "analysis/ADRD_agg_lagged_trimmed_1_99.fst"))


## Formulas for outcome models (parametric and semiparametric thin-plate spline)

formula_expos_only <- as.formula(paste(outcome_name, "~", paste(c("w", strata_vars), collapse = "+", sep = "")))
formula_expos_only_smooth <- as.formula(paste(outcome_name, "~", paste(c("s(w, bs = 'ts')", strata_vars), collapse = "+", sep = "")))