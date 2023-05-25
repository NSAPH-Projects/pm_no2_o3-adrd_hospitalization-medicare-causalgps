## User should update directory paths for data, code, and results ##

dir_data <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/data/"
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/code/"
dir_results <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/results/"

## Classify variables in dataset

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


## Formulas for outcome models (parametric and semiparametric thin-plate spline)

formula_expos_only <- as.formula(paste("Y ~", paste(c("w", strata_vars), collapse = "+", sep = "")))
formula_expos_only_smooth <- as.formula(paste("Y ~", paste(c("s(w, bs = 'ts')", strata_vars), collapse = "+", sep = "")))
formula_expos_only_smooth_cr <- as.formula(paste("Y ~", paste(c("s(w, bs = 'cr', k = 4, fx = TRUE)", strata_vars), collapse = "+", sep = "")))
