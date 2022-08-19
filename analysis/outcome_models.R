library(mgcv)

##### Naive (associational) Poisson regression #####

bam_naive <- bam(Y_subset ~ w_subset + no2 + ozone_summer +
                   any_dual + ADRD_age + sexM + race_cat +
                   summer_tmmx + summer_rmax + region + ADRD_year +
                   mean_bmi + smoke_rate + hispanic + pct_blk +
                   medhouseholdincome + medianhousevalue + PIR + poverty +
                   education + popdensity + pct_owner_occ,
                 data = data_subset,
                 offset = log(person_years),
                 family = poisson(link = "log"),
                 samfrac = 0.05,
                 chunk.size = 5000,
                 control = gam.control(trace = TRUE))
summary(bam_naive)


##### Poisson regression adjusting for GPS #####

data_with_gps <- copy(data_subset)
data_with_gps$gps <- estimated_gps$gps

bam_gps_adjusted <- bam(Y_subset ~ w_subset + no2 + ozone_summer +
                   any_dual + ADRD_age + sexM + race_cat +
                   summer_tmmx + summer_rmax + region + ADRD_year +
                   mean_bmi + smoke_rate + hispanic + pct_blk +
                   medhouseholdincome + medianhousevalue + PIR + poverty +
                   education + popdensity + pct_owner_occ + gps,
                 data = data_with_gps,
                 offset = log(person_years),
                 family = poisson(link = "log"),
                 samfrac = 0.05,
                 chunk.size = 5000,
                 control = gam.control(trace = TRUE))
summary(bam_gps_adjusted)


##### (Poisson) Parametric outcome models #####

# method 1: gam package
# To do: rd rest of Kevin's code (https://github.com/kevjosey/erc-strata/blob/master/R/match_estimate.R)
# to see if there's anything extra I should do to account for matching

# GPS matching
bam_exposure_only_matched <- bam(Y ~ w + any_dual + ADRD_age + sexM + race_cat,
                                 data = matched_data,
                                 offset = log(person_years),
                                 family = poisson(link = "log"),
                                 weights = counter,
                                 samfrac = 0.05,
                                 chunk.size = 5000,
                                 control = gam.control(trace = TRUE))
summary(bam_exposure_only_matched)

bam_exposures_controlled_matched <- bam(Y ~ w + no2 + ozone_summer +
                                          any_dual + ADRD_age + sexM + race_cat,
                                        data = matched_data,
                                        offset = log(person_years),
                                        family = poisson(link = "log"),
                                        weights = counter,
                                        samfrac = 0.05,
                                        chunk.size = 5000,
                                        control = gam.control(trace = TRUE))
summary(bam_exposures_controlled_matched)

bam_doubly_robust_matched <- bam(Y ~ w + no2 + ozone_summer +
                                   any_dual + ADRD_age + sexM + race_cat +
                                   summer_tmmx + summer_rmax + region + ADRD_year +
                                   mean_bmi + smoke_rate + hispanic + pct_blk +
                                   medhouseholdincome + medianhousevalue + PIR + poverty +
                                   education + popdensity + pct_owner_occ,
                                 data = matched_data,
                                 offset = log(person_years),
                                 family = poisson(link = "log"),
                                 weights = counter,
                                 samfrac = 0.05,
                                 chunk.size = 5000,
                                 control = gam.control(trace = TRUE))
summary(bam_doubly_robust_matched)


### GPS weighting
bam_exposure_only_weighted <- bam(Y ~ w + any_dual + ADRD_age + sexM + race_cat,
                                  data = weighted_data,
                                  offset = log(person_years),
                                  family = poisson(link = "log"),
                                  weights = ipw,
                                  samfrac = 0.05,
                                  chunk.size = 5000,
                                  control = gam.control(trace = TRUE))
summary(bam_exposure_only_weighted)

bam_exposures_controlled_weighted <- bam(Y ~ w + no2 + ozone_summer +
                                           any_dual + ADRD_age + sexM + race_cat,
                                         data = weighted_data,
                                         offset = log(person_years),
                                         family = poisson(link = "log"),
                                         weights = ipw,
                                         samfrac = 0.05,
                                         chunk.size = 5000,
                                         control = gam.control(trace = TRUE))
summary(bam_exposures_controlled_weighted)

bam_doubly_robust_weighted <- bam(Y ~ w + no2 + ozone_summer +
                                    any_dual + ADRD_age + sexM + race_cat +
                                    summer_tmmx + summer_rmax + region + ADRD_year +
                                    mean_bmi + smoke_rate + hispanic + pct_blk +
                                    medhouseholdincome + medianhousevalue + PIR + poverty +
                                    education + popdensity + pct_owner_occ,
                                  data = weighted_data,
                                  offset = log(person_years),
                                  family = poisson(link = "log"),
                                  weights = ipw,
                                  samfrac = 0.05,
                                  chunk.size = 5000,
                                  control = gam.control(trace = TRUE))
summary(bam_doubly_robust_weighted)

### Capped GPS weighting
bam_exposure_only_weighted <- bam(Y ~ w + any_dual + ADRD_age + sexM + race_cat,
                                  data = capped_weighted_data,
                                  offset = log(person_years),
                                  family = poisson(link = "log"),
                                  weights = ipw,
                                  samfrac = 0.05,
                                  chunk.size = 5000,
                                  control = gam.control(trace = TRUE))
summary(bam_exposure_only_weighted)

bam_exposures_controlled_weighted <- bam(Y ~ w + no2 + ozone_summer +
                                           any_dual + ADRD_age + sexM + race_cat,
                                         data = capped_weighted_data,
                                         offset = log(person_years),
                                         family = poisson(link = "log"),
                                         weights = ipw,
                                         samfrac = 0.05,
                                         chunk.size = 5000,
                                         control = gam.control(trace = TRUE))
summary(bam_exposures_controlled_weighted)

bam_doubly_robust_weighted <- bam(Y ~ w + no2 + ozone_summer +
                                    any_dual + ADRD_age + sexM + race_cat +
                                    summer_tmmx + summer_rmax + region + ADRD_year +
                                    mean_bmi + smoke_rate + hispanic + pct_blk +
                                    medhouseholdincome + medianhousevalue + PIR + poverty +
                                    education + popdensity + pct_owner_occ,
                                  data = capped_weighted_data,
                                  offset = log(person_years),
                                  family = poisson(link = "log"),
                                  weights = ipw,
                                  samfrac = 0.05,
                                  chunk.size = 5000,
                                  control = gam.control(trace = TRUE))
summary(bam_doubly_robust_weighted)


# method 2: gnm package
# To do: fix error: only first element used
# rd rest of Xiao's code (https://github.com/wxwx1993/National_Causal/blob/master/statistical_models.R)
gnm_doubly_robust <- gnm(Y ~ w + no2 + ozone_summer +
                           any_dual + ADRD_age + sexM + race_cat +
                           summer_tmmx + summer_rmax + region + ADRD_year +
                           mean_bmi + smoke_rate + hispanic + pct_blk +
                           medhouseholdincome + medianhousevalue + PIR + poverty +
                           education + popdensity + pct_owner_occ,
                         eliminate = (as.factor(sexM):as.factor(race_cat):as.factor(any_dual):ADRD_age),
                         data = matched_data,
                         offset = log(person_years),
                         family = poisson(link = "log"),
                         weights = counter)
summary(gnm_doubly_robust)

# method 3: CausalGPS package
# To do: wrong right now cuz offset should be constant not various in the regression; to do - fix

outcome <- estimate_pmetric_erf(formula = Y ~ . -row_index -gps -counter, # w or .?
                                family = poisson, # poisson(link = "log")
                                data = matched_data, # To Do: check if this must be a data.frame
                                ci_appr = "matching")
summary(outcome)




##### Print results from (Poisson) parametric outcome models #####

iqr <- IQR(ADRD_agg_lagged$pm25)
iqr
coef_bam_doubly_robust <- bam_doubly_robust$coefficients[2]
coef_bam_exposures_controlled <- bam_exposures_controlled$coefficients[2]
coef_bam_exposure_only <- bam_exposure_only$coefficients[2]
coef_bam_doubly_robust
coef_bam_exposures_controlled
coef_bam_exposure_only
cat("Hazard ratio per 1 IQR increase in PM2.5:", exp(coef_bam_doubly_robust*iqr))
cat("Hazard ratio per 1 IQR increase in PM2.5:", exp(coef_bam_exposures_controlled*iqr))
cat("Hazard ratio per 1 IQR increase in PM2.5:", exp(coef_bam_exposure_only*iqr))


##### Non-parametric model on matched and weighted data ##### 

# model rate log(Y/offset) instead of Y, to incorporate offset
matched_erf <- estimate_npmetric_erf(as.double(log(matched_data$Y / matched_data$person_years + 0.001)),
                                     as.double(matched_data$w),
                                     matched_data$counter,
                                     bw_seq=seq(0.1, 10, length.out = 15),
                                     w_vals = seq(0, range(matched_data$w)[2], length.out = 100),
                                     nthread = 16)
matched_erf_plot <- plot(matched_erf) # plot log_rate
matched_erf_plot$erf <- exp(matched_erf_plot$erf)
plot(matched_erf) # plot rate

weighted_erf <- estimate_npmetric_erf(as.double(log(weighted_data$Y / weighted_data$person_years + 0.001)),
                                      as.double(weighted_data$w),
                                      round(weighted_data$ipw * 50),
                                      bw_seq=seq(0.1, 10, length.out = 15),
                                      w_vals = seq(0, range(weighted_data$w)[2], length.out = 100),
                                      nthread = 16)
weighted_erf_plot <- plot(weighted_erf) # plot log_rate
weighted_erf_plot$erf <- exp(weighted_erf_plot$erf)
plot(matched_erf) # plot rate
#summary(weighted_erf)


##### Semi-parametric model on matched and weighted data ##### 

# to do: fix error
# to do: try Y~s(w,df=3)
semipmetric_exposure_only <- estimate_semipmetric_erf(log(Y/person_years + 0.001) ~ w + any_dual + ADRD_age + sexM + race_cat,
                                                      data = matched_data,
                                                      family = poisson(link = "log"),
                                                      ci_appr = "matching")
matched_semi_erf_plot <- plot(semipmetric_exposure_only) # plot log_rate
matched_semi_erf_plot$erf <- exp(matched_semi_erf_plot$erf)
plot(matched_semi_erf_plot)



#### To Do: use transformed variables created in explore_covariate_distributions.R ####


##### Old Code #####

# Outcome models

model_pm25 <- bam(n_ADRDhosp ~ offset(log(n_persons * n_years)) +
                    pm25 + no2 + ozone_summer + tmmx + rmax +
                    factor(ffs_entry_year) + factor(n_years) + factor(any_dual) +
                    ADRD_age + factor(sexM) + factor(race_cat) +
                    mean_bmi + smoke_rate + hispanic + pct_blk +
                    PIR + poverty +
                    education + popdensity + pct_owner_occ,
                  data = ADRD_agg_lagged,
                  weights = ADRD_agg_lagged[, ipw_pm25],
                  family = poisson, 
                  samfrac = 0.05, chunk.size = 5000,
                  control = gam.control(trace = TRUE))

model_no2 <- bam(n_ADRDhosp ~ offset(log(n_persons * n_years)) +
                   pm25 + no2 + ozone_summer + tmmx + rmax +
                   factor(ffs_entry_year) + factor(n_years) + factor(any_dual) +
                   ADRD_age + factor(sexM) + factor(race_cat) +
                   mean_bmi + smoke_rate + hispanic + pct_blk +
                   PIR + poverty +
                   education + popdensity + pct_owner_occ,
                 data = ADRD_agg_lagged,
                 weights = ADRD_agg_lagged[, ipw_no2],
                 family = poisson, 
                 samfrac = 0.05, chunk.size = 5000,
                 control = gam.control(trace = TRUE))

model_ozone <- bam(n_ADRDhosp ~ offset(log(n_persons * n_years)) +
                     pm25 + no2 + ozone_summer + tmmx + rmax +
                     factor(ffs_entry_year) + factor(n_years) + factor(any_dual) +
                     ADRD_age + factor(sexM) + factor(race_cat) +
                     mean_bmi + smoke_rate + hispanic + pct_blk +
                     PIR + poverty +
                     education + popdensity + pct_owner_occ,
                   data = ADRD_agg_lagged,
                   weights = ADRD_agg_lagged[, ipw_ozone],
                   family = poisson, 
                   samfrac = 0.05, chunk.size = 5000,
                   control = gam.control(trace = TRUE))

# Results

cat("IQR of PM2.5 (micrograms/m^3):", IQR(ADRD_agg_lagged$pm25))
cat("IQR of NO2 (ppb):", IQR(ADRD_agg_lagged$no2))
cat("IQR of summer ozone (ppb):", IQR(ADRD_agg_lagged$ozone_summer))

cat("Hazard ratio per 1 IQR increase in PM2.5:", exp(model_pm25$coefficients[2]*IQR(ADRD_agg_lagged$pm25)))
cat("Hazard ratio per 1 IQR increase in NO2:", exp(model_no2$coefficients[2]*IQR(ADRD_agg_lagged$no2)))
cat("Hazard ratio per 1 IQR increase in summer ozone:", exp(model_ozone$coefficients[2]*IQR(ADRD_agg_lagged$ozone_summer)))

print("To compare with Shi et al 2021 preprint")
cat("Hazard ratio per 3.2 micrograms/m^3 increase in PM2.5:", exp(model_pm25$coefficients[2]*3.2))
cat("Hazard ratio per 11.6 ppb increase in NO2:", exp(model_no2$coefficients[2]*11.6))
cat("Hazard ratio per 5.3 ppb increase in summer ozone:", exp(model_ozone$coefficients[2]*5.3))

print("To compare with Shi et al 2020")
cat("Hazard ratio per 5 micrograms/m^3 increase in PM2.5:", exp(model_pm25$coefficients[2]*5))


# model.matrix(y~-1+as.factor(year)+age,data)
# c <- as.matrix(c) # don't use this
# c <- as.data.frame(c)

