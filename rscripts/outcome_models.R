#TODO: Convert it to project's convention.

library(mgcv)

# Formulas for Poisson regression
formula_all_covars <- as.formula(paste("Y ~", paste(c("w", indiv_var_names, other_expos_names, zip_var_names), collapse = "+", sep = "")))
formula_all_expos <- as.formula(paste("Y ~", paste(c("w", indiv_var_names, other_expos_names), collapse = "+", sep = "")))
formula_expos_only <- as.formula(paste("Y ~", paste(c("w", indiv_var_names), collapse = "+", sep = "")))

# Formulas for Poisson regression, with thin-plate spline for exposure
formula_all_covars_smooth <- as.formula(paste("Y ~", paste(c("s(w, bs = 'ts')", indiv_var_names, other_expos_names, zip_var_names), collapse = "+", sep = "")))
formula_all_expos_smooth <- as.formula(paste("Y ~", paste(c("s(w, bs = 'ts')", indiv_var_names, other_expos_names), collapse = "+", sep = "")))
formula_expos_only_smooth <- as.formula(paste("Y ~", paste(c("s(w, bs = 'ts')", indiv_var_names), collapse = "+", sep = "")))


##### Naive (associational) Poisson regression #####

bam_naive_expos_only <- bam(formula_expos_only,
                            data = all_data,
                            offset = log(person_years),
                            family = poisson(link = "log"),
                            samfrac = 0.05,
                            chunk.size = 5000,
                            control = gam.control(trace = TRUE))
summary(bam_naive_expos_only)
saveRDS(summary(bam_naive_expos_only), file = paste0(dir_proj, "results/parametric_results/bam_naive_exposure_only_", n_rows, "rows_", modifications, ".rds"))

bam_naive_all_expos <- bam(formula_all_expos,
                 data = all_data,
                 offset = log(person_years),
                 family = poisson(link = "log"),
                 samfrac = 0.05,
                 chunk.size = 5000,
                 control = gam.control(trace = TRUE))
summary(bam_naive_all_expos)
saveRDS(summary(bam_naive_all_expos), file = paste0(dir_proj, "results/parametric_results/bam_naive_all_exposures_", n_rows, "rows_", modifications, ".rds"))

bam_naive_all_covars <- bam(formula_all_covars,
                            data = all_data,
                            offset = log(person_years),
                            family = poisson(link = "log"),
                            samfrac = 0.05,
                            chunk.size = 5000,
                            control = gam.control(trace = TRUE))
summary(bam_naive_all_covars)
saveRDS(summary(bam_naive_all_covars), file = paste0(dir_proj, "results/parametric_results/bam_naive_all_covariates_", n_rows, "rows_", modifications, ".rds"))



##### Smoothed naive (associational) Poisson regression #####

bam_smooth_naive_expos_only <- bam(formula_expos_only_smooth,
                            data = all_data,
                            offset = log(person_years),
                            family = poisson(link = "log"),
                            samfrac = 0.05,
                            chunk.size = 5000,
                            control = gam.control(trace = TRUE),
                            nthreads = n_cores - 1)
png(paste0(dir_proj, "results/semiparametric_results/ERFs/bam_smooth_naive_exposure_only_", n_rows, "rows_", modifications, ".png"))
plot(bam_smooth_naive_expos_only, main = paste0("Naive Poisson, exposure only (", exposure_name, ")"))
dev.off()

bam_smooth_naive_all_expos <- bam(formula_all_expos_smooth,
                           data = all_data,
                           offset = log(person_years),
                           family = poisson(link = "log"),
                           samfrac = 0.05,
                           chunk.size = 5000,
                           control = gam.control(trace = TRUE),
                           nthreads = n_cores - 1)
png(paste0(dir_proj, "results/semiparametric_results/ERFs/bam_smooth_naive_all_exposures_", n_rows, "rows_", modifications, ".png"))
plot(bam_smooth_naive_all_expos, main = paste("Naive Poisson,", exposure_name, "as exposure, all exposures"))
dev.off()

bam_smooth_naive_all_covars <- bam(formula_all_covars_smooth,
                            data = all_data,
                            offset = log(person_years),
                            family = poisson(link = "log"),
                            samfrac = 0.05,
                            chunk.size = 5000,
                            control = gam.control(trace = TRUE),
                            nthreads = n_cores - 1)
png(paste0(dir_proj, "results/semiparametric_results/ERFs/bam_smooth_naive_all_covariates_", n_rows, "rows_", modifications, ".png"))
plot(bam_smooth_naive_all_covars, main = paste("Naive Poisson,", exposure_name, "as exposure, all covariates"))
dev.off()


##### Poisson regression matching on GPS #####

# method 1: gam package
# To do: rd rest of Kevin's code (https://github.com/kevjosey/erc-strata/blob/master/R/match_estimate.R)
# to see if there's anything extra I should do to account for matching

# GPS matching
bam_exposure_only_matched <- bam(formula_expos_only,
                                 data = matched_data,
                                 offset = log(person_years),
                                 family = poisson(link = "log"),
                                 weights = counter_weight,
                                 samfrac = 0.05,
                                 chunk.size = 5000,
                                 control = gam.control(trace = TRUE))
summary(bam_exposure_only_matched)
saveRDS(summary(bam_exposure_only_matched), file = paste0(dir_proj, "results/parametric_results/bam_matched_exposure_only_", n_rows, "rows_", modifications, ".rds"))

bam_exposures_controlled_matched <- bam(formula_all_expos,
                                        data = matched_data,
                                        offset = log(person_years),
                                        family = poisson(link = "log"),
                                        weights = counter_weight,
                                        samfrac = 0.05,
                                        chunk.size = 5000,
                                        control = gam.control(trace = TRUE))
summary(bam_exposures_controlled_matched)
saveRDS(summary(bam_exposures_controlled_matched), file = paste0(dir_proj, "results/parametric_results/bam_matched_all_exposures_", n_rows, "rows_", modifications, ".rds"))

bam_all_covariates_matched <- bam(formula_all_covars,
                                 data = matched_data,
                                 offset = log(person_years),
                                 family = poisson(link = "log"),
                                 weights = counter_weight,
                                 samfrac = 0.05,
                                 chunk.size = 5000,
                                 control = gam.control(trace = TRUE))
summary(bam_all_covariates_matched)
saveRDS(summary(bam_all_covariates_matched), file = paste0(dir_proj, "results/parametric_results/bam_matched_all_covariates_", n_rows, "rows_", modifications, ".rds"))

##### Smoothed Poisson regression matching on GPS

bam_smooth_exposure_only_matched <- bam(formula_expos_only_smooth,
                                                data = matched_data,
                                                offset = log(person_years),
                                                family = poisson(link = "log"),
                                                weights = counter_weight,
                                                samfrac = 0.05,
                                                chunk.size = 5000,
                                                control = gam.control(trace = TRUE),
                                                nthreads = n_cores - 1)
png(paste0(dir_proj, "results/semiparametric_results/ERFs/bam_smooth_exposure_only_matched_", n_rows, "rows_", modifications, ".png"))
plot(bam_smooth_exposure_only_matched, main = paste0("GPS-Matched, Poisson regression,\nexposure only (", exposure_name, ")"))
dev.off()

bam_smooth_exposures_controlled_matched <- bam(formula_all_expos_smooth,
                                                       data = matched_data,
                                                       offset = log(person_years),
                                                       family = poisson(link = "log"),
                                                       weights = counter_weight,
                                                       samfrac = 0.05,
                                                       chunk.size = 5000,
                                                       control = gam.control(trace = TRUE),
                                                       nthreads = n_cores - 1)
png(paste0(dir_proj, "results/semiparametric_results/ERFs/bam_smooth_all_expos_matched_", n_rows, "rows_", modifications, ".png"))
plot(bam_smooth_exposures_controlled_matched, main = paste0("GPS-Matched, Poisson regression\n", exposure_name, " as exposure, other exposures controlled"))
dev.off()

bam_smooth_all_covariates_matched <- bam(formula_all_covars_smooth,
                                                 data = matched_data,
                                                 offset = log(person_years),
                                                 family = poisson(link = "log"),
                                                 weights = counter_weight,
                                                 samfrac = 0.05,
                                                 chunk.size = 5000,
                                                 control = gam.control(trace = TRUE),
                                                 nthreads = n_cores - 1)
png(paste0(dir_proj, "results/semiparametric_results/ERFs/bam_smooth_all_covars_matched_", n_rows, "rows_", modifications, ".png"))
plot(bam_smooth_all_covariates_matched, main = paste0("GPS-Matched, Poisson regression\n", exposure_name, " as exposure, all covariates controlled"))
dev.off()


##### Poisson regression weighting (CAPPED) on GPS #####
bam_exposure_only_capped_weighted <- bam(formula_expos_only,
                                  data = capped_weighted_data,
                                  offset = log(person_years),
                                  family = poisson(link = "log"),
                                  weights = counter_weight,
                                  samfrac = 0.05,
                                  chunk.size = 5000,
                                  control = gam.control(trace = TRUE))
summary(bam_exposure_only_capped_weighted)
saveRDS(summary(bam_exposure_only_capped_weighted), file = paste0(dir_proj, "results/parametric_results/bam_capped_weighted_exposure_only_", n_rows, "rows_", modifications, ".rds"))

bam_exposures_controlled_capped_weighted <- bam(formula_all_expos,
                                         data = capped_weighted_data,
                                         offset = log(person_years),
                                         family = poisson(link = "log"),
                                         weights = counter_weight,
                                         samfrac = 0.05,
                                         chunk.size = 5000,
                                         control = gam.control(trace = TRUE))
summary(bam_exposures_controlled_capped_weighted)
saveRDS(summary(bam_exposures_controlled_capped_weighted), file = paste0(dir_proj, "results/parametric_results/bam_capped_weighted_all_exposures_", n_rows, "rows_", modifications, ".rds"))

bam_all_covariates_capped_weighted <- bam(formula_all_covars,
                                  data = capped_weighted_data,
                                  offset = log(person_years),
                                  family = poisson(link = "log"),
                                  weights = counter_weight,
                                  samfrac = 0.05,
                                  chunk.size = 5000,
                                  control = gam.control(trace = TRUE))
summary(bam_all_covariates_capped_weighted)
saveRDS(summary(bam_all_covariates_capped_weighted), file = paste0(dir_proj, "results/parametric_results/bam_weighted_all_covariates_", n_rows, "rows_", modifications, ".rds"))


##### Smoothed Poisson regression weighting (CAPPED) on GPS

bam_smooth_exposure_only_capped_weighted <- bam(formula_expos_only_smooth,
                                   data = capped_weighted_data,
                                   offset = log(person_years),
                                   family = poisson(link = "log"),
                                   weights = counter_weight,
                                   samfrac = 0.05,
                                   chunk.size = 5000,
                                   control = gam.control(trace = TRUE),
                                   nthreads = n_cores - 1)
png(paste0(dir_proj, "results/semiparametric_results/ERFs/bam_smooth_exposure_only_capped_weighted_", n_rows, "rows_", modifications, ".png"))
plot(bam_smooth_exposure_only_capped_weighted, main = paste0("GPS-Weighted, Capped at 10, Poisson regression,\nexposure only (", exposure_name, ")"))
dev.off()
saveRDS(bam_smooth_exposure_only_capped_weighted, file = paste0(dir_proj, "results/semiparametric_results/spline_objects/bam_smooth_exposure_only_capped_weighted_", n_rows, "rows_", modifications, ".rds"))

bam_smooth_exposures_controlled_capped_weighted <- bam(formula_all_expos_smooth,
                                                data = capped_weighted_data,
                                                offset = log(person_years),
                                                family = poisson(link = "log"),
                                                weights = counter_weight,
                                                samfrac = 0.05,
                                                chunk.size = 5000,
                                                control = gam.control(trace = TRUE),
                                                nthreads = n_cores - 1)
png(paste0(dir_proj, "results/semiparametric_results/ERFs/bam_smooth_all_expos_capped_weighted_", n_rows, "rows_", modifications, ".png"))
plot(bam_smooth_exposures_controlled_capped_weighted, main = paste0("GPS-Weighted, Capped at 10, Poisson regression\n", exposure_name, " as exposure, other exposures controlled"))
dev.off()
saveRDS(bam_smooth_exposures_controlled_capped_weighted, file = paste0(dir_proj, "results/semiparametric_results/spline_objects/bam_smooth_all_expos_capped_weighted_", n_rows, "rows_", modifications, ".rds"))

bam_smooth_all_covariates_capped_weighted <- bam(formula_all_covars_smooth,
                                                       data = capped_weighted_data,
                                                       offset = log(person_years),
                                                       family = poisson(link = "log"),
                                                       weights = counter_weight,
                                                       samfrac = 0.05,
                                                       chunk.size = 5000,
                                                       control = gam.control(trace = TRUE),
                                                       nthreads = n_cores - 1)
png(paste0(dir_proj, "results/semiparametric_results/ERFs/bam_smooth_all_covars_capped_weighted_", n_rows, "rows_", modifications, ".png"))
plot(bam_smooth_all_covariates_capped_weighted, main = paste0("GPS-Weighted, Capped at 10, Poisson regression\n", exposure_name, " as exposure, all covariates controlled"))
dev.off()
saveRDS(bam_smooth_all_covariates_capped_weighted, file = paste0(dir_proj, "results/semiparametric_results/spline_objects/bam_smooth_all_covars_capped_weighted_", n_rows, "rows_", modifications, ".rds"))



##### Old code #####

# cat("IQR of PM2.5 (micrograms/m^3):", IQR(ADRD_agg_lagged$pm25))
# cat("IQR of NO2 (ppb):", IQR(ADRD_agg_lagged$no2))
# cat("IQR of summer ozone (ppb):", IQR(ADRD_agg_lagged$ozone_summer))

# cat("Hazard ratio per 1 IQR increase in PM2.5:", exp(model_pm25$coefficients[2]*IQR(ADRD_agg_lagged$pm25)))
# cat("Hazard ratio per 1 IQR increase in NO2:", exp(model_no2$coefficients[2]*IQR(ADRD_agg_lagged$no2)))
# cat("Hazard ratio per 1 IQR increase in summer ozone:", exp(model_ozone$coefficients[2]*IQR(ADRD_agg_lagged$ozone_summer)))

# print("To compare with Shi et al 2021 preprint")
# cat("Hazard ratio per 3.2 micrograms/m^3 increase in PM2.5:", exp(model_pm25$coefficients[2]*3.2))
# cat("Hazard ratio per 11.6 ppb increase in NO2:", exp(model_no2$coefficients[2]*11.6))
# cat("Hazard ratio per 5.3 ppb increase in summer ozone:", exp(model_ozone$coefficients[2]*5.3))

# print("To compare with Shi et al 2020")
# cat("Hazard ratio per 5 micrograms/m^3 increase in PM2.5:", exp(model_pm25$coefficients[2]*5))



