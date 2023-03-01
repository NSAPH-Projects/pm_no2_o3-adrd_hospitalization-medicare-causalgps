offset_var_names <- c("n_persons", "n_years")
zip_expos_names <- c("pm25", "no2", "ozone_summer")
zip_quant_var_names <- c("mean_bmi", "smoke_rate", "hispanic", "prop_blk",
                         "PIR", "poverty", "education", "popdensity", "prop_owner_occ",
                         "summer_tmmx", "summer_rmax")
zip_unordered_cat_var_names <- c("region")
indiv_quant_var_names <- NULL
indiv_unordered_cat_var_names <- c("sex", "race", "dual", "age_grp", "year") # to do: age_group is ordered
strata_vars <- c("year", "sex", "race", "dual", "age_grp") # to do: consider if year should be here
zip_var_names <- c(zip_quant_var_names, zip_unordered_cat_var_names)
indiv_var_names <- c(indiv_unordered_cat_var_names, indiv_quant_var_names) # note: for now, using ADRD_age as a quantitative variable (not binned)



