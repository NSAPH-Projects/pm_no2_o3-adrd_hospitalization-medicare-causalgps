### NOTE: User should update directory paths for data, code, and results ###

# directories for code and results to be saved (separate from data directory)
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/git/code/"
dir_results <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/git/results/"

# directory for cleaned, processed, and aggregated data (PRIVATE CMS DATA, SHOULD NOT BE EXPORTED)
dir_data <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/data/"

# directories for private raw data (PRIVATE CMS DATA, SHOULD NOT BE EXPORTED)
dir_denominator <- "/n/dominici_nsaph_l3/Lab/projects/analytic/denom_by_year/"
dir_hosp <- "/n/dominici_nsaph_l3/Lab/projects/analytic/adrd_hospitalization/"

# directories for public raw data
dir_pm25 <- "/nfs/home/D/dam9096/shared_space/ci3_exposure/pm25/whole_us/annual/zipcode/qd_predictions_ensemble/ywei_aggregation/"
dir_pm25_comp <- "/nfs/home/D/dam9096/shared_space/ci3_exposure/pm25_components/"
dir_no2 <- "/nfs/home/D/dam9096/shared_space/ci3_exposure/no2/whole_us/annual/zipcode/qd_predictions_ensemble/ywei_aggregations/"
dir_ozone <- "/nfs/home/D/dam9096/shared_space/ci3_exposure/ozone/whole_us/seasonal/zipcode/"
dir_temp <- "/nfs/home/D/dam9096/shared_space/ci3_confounders/data_for_analysis/prepped_temperature/annual/"
dir_confounders <- "/nfs/home/D/dam9096/shared_space/ci3_analysis/National_Causal/National_Causal/2016_temp/aggregate_data.RData"


### Variables and formulas for this analysis (after the data have been cleaned, processed, and aggregated) ###

# Classify variables in dataset

offset_var_names <- c("n_persons", "n_years")
zip_expos_names <- c("pm25", "no2", "ozone_summer")
zip_quant_var_names <- c("mean_bmi", "smoke_rate", "hispanic", "prop_blk",
                         "PIR", "poverty", "education", "popdensity", "prop_owner_occ",
                         "summer_tmmx", "summer_rmax")
zip_unordered_cat_var_names <- c("region")
indiv_quant_var_names <- NULL
indiv_unordered_cat_var_names <- c("sex", "race", "dual", "age_grp", "year") # note that age_grp here is an unordered categorical variable
strata_vars <- c("year", "sex", "race", "dual", "age_grp") # include year to account for possible confounding, such as billing incentives changing year to year
zip_var_names <- c(zip_quant_var_names, zip_unordered_cat_var_names)
indiv_var_names <- c(indiv_unordered_cat_var_names, indiv_quant_var_names)

# outcome variable for this analysis
outcome_name <- "n_hosp"


# Formulas for outcome models (parametric and semiparametric thin-plate spline)

formula_expos_only <- as.formula(paste("Y ~", paste(c("w", strata_vars), collapse = "+", sep = "")))
formula_expos_only_smooth <- as.formula(paste("Y ~", paste(c("s(w, bs = 'ts')", strata_vars), collapse = "+", sep = "")))
formula_expos_only_smooth_cr <- as.formula(paste("Y ~", paste(c("s(w, bs = 'cr', k = 4, fx = TRUE)", strata_vars), collapse = "+", sep = "")))
