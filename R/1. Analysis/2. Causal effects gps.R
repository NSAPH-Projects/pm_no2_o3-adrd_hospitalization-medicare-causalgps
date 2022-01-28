###############################################################################
# Project: Causal effect of air pollution on first ADRD hospitalization       #
# Code: estimate the causal effect of air pollution using GPS weighting       #
# Input: "ADRD_aggregated_entry_to_followup_merged.fst"                       #                 
# Author: Daniel Mork                                                         #
# Date: 2021-12-29                                                            #
###############################################################################
rm(list = ls())
gc()
############################# 0. Setup ########################################
library(NSAPHutils)
library(data.table)
library(fst)
library(zipcode)
library(CausalGPS)
library(mgcv)
setDTthreads(threads = 16)

dir_data <- "/nfs/home/D/dam9096/shared_space/ci3_analysis/dmork/Data/Causal_ADRD/"
dir_gps <- "/nfs/home/D/dam9096/shared_space/ci3_analysis/dmork/Data/Causal_ADRD/GPS_Data/"

############################# 1. Load data #####################################
ADRDdat <- read_fst(paste0(dir_data, "ADRD_aggregated_entry_to_followup_merged.fst"), as.data.table = TRUE)
names(ADRDdat)
# [1] "followupZip"  "followupYear" "sex"          "race"         "dual"        
# [6] "entryYear"    "age"          "nADRDEvents"  "nDeaths"      "N"           
# [11] "nYears"       "followupAge"  "avg_age"      "med_age"      "bmi"         
# [16] "smoke"        "hispanic"     "pct_blk"      "income"       "value"       
# [21] "poverty"      "education"    "popdensity"   "ownerocc"     "summer_tmmx" 
# [26] "summer_rmax"  "winter_tmmx"  "winter_rmax"  "death_prob"   "pm25"        
# [31] "no2"          "ozone"        "nox"          "city"         "state"       
# [36] "latitude"     "longitude"    "region"      
ADRDdat <- ADRDdat[complete.cases(ADRDdat)]
ADRDdat[sex == 2, sex := 0] # binary, 1 = male
ADRDdat$followupYear <- as.factor(ADRDdat$followupYear)


######################### 2. Estimate GPS weights ##############################
ADRDgps <- estimate_gps(ADRDdat$nADRDEvents, ADRDdat$pm25, 
                        ADRDdat[, .(sex, race, dual, followupYear,
                                    bmi, smoke, hispanic, pct_blk,
                                    income, value, poverty, education,
                                    popdensity, ownerocc,
                                    summer_tmmx, winter_tmmx, 
                                    summer_rmax, winter_rmax, region)],
                        pred_model = "sl",
                        gps_model = "non-parametric",
                        internal_use = TRUE,
                        params = list(xgb_rounds = 5), # increase later...
                        nthread = 16,                                
                        sl_lib = c("m_xgboost"))
write_fst(ADRDgps, paste0(dir_gps, "ADRDgps.fst"))

######################### 3. Check covariate balance ###########################
bal <- CausalGPS::absolute_weighted_corr_fun(
  ADRDdat$pm25,
  ADRDgps[[1]]$gps,
  ADRDdat[, .(sex, race, dual, followupYear,
              bmi, smoke, hispanic, pct_blk,
              income, value, poverty, education,
              popdensity, ownerocc,
              summer_tmmx, winter_tmmx, 
              summer_rmax, winter_rmax,
              latitude, longitude)])





######################### 4. Calculate IPW #####################################
ADRD_ipw <- dnorm(ADRDdat$pm25, mean(ADRDdat$pm25), sd(ADRDdat$pm25)) / 
  ADRDgps[[1]]$gps


######################### 5. Linear GLM fit ####################################
m1 <- glm(nADRDEvents ~ 
            offset(log(I(nYears * N))) + # offset total person-time at-risk
            pm25,
          data = ADRDdat, weights = ADRD_ipw, 
          family = poisson(link = "log")) # print diagnostic output
summary(m1)
exp(sort(coef(m1))) # sorted effects




######################### 5. Smooth GAM fit ####################################
m1 <- bam(N ~ offset(log(size)) +
            s(pm25, bs = "cr", k = 3) +
            s(ox, bs = "cr", k = 3) +
            region + as.factor(year) +
            poly(pct_dual, 2) + 
            pct_male + pct_oth_race + pct_blk + pct_his +
            poly(age, 2) + poly(bmi, 2) + 
            smoke + income + value +
            poverty + ed + dens + ownocc, 
          data = ADRD1, weights = ADRD_pm25_gps[[1]]$gps,
          family = poisson(link = "log"),
          discrete = TRUE,
          nthreads = 8)
summary(m1)
plot(m1, select = 1, scale = 0, trans = exp)
plot(m1, select = 2, scale = 0, trans = exp)

pm25_plot <- plot(m1, select = 1, xlim = quantile(ADRD1$pm25, c(0.005, 0.995), na.rm = T))
plot(m1, select = 1, trans = exp, scale = 0, shift = -(pm25_plot[[1]]$fit[1]),
     xlim = quantile(ADRD1$pm25, c(0.005, 0.995), na.rm = T))
abline(v = quantile(ADRD1$pm25, c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = T))
abline(h = exp(pm25_plot[[1]]$fit[sapply(c(0.1, 0.25, 0.5, 0.75, 0.9), function(q) {
  which.min(abs(pm25_plot[[1]]$x - quantile(ADRD1$pm25, q, na.rm = T)))
})] - (pm25_plot[[1]]$fit[1])))


ox_plot <- plot(m1, select = 2, xlim = quantile(ADRD1$ox, c(0.005, 0.995), na.rm = T))
plot(m1, select = 2, trans = exp, scale = 0, shift = -(ox_plot[[2]]$fit[1]),
     xlim = quantile(ADRD1$ox, c(0.005, 0.995), na.rm = T))
abline(v = quantile(ADRD1$ox, c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = T))
abline(h = exp(ox_plot[[2]]$fit[sapply(c(0.1, 0.25, 0.5, 0.75, 0.9), function(q) {
  which.min(abs(ox_plot[[2]]$x - quantile(ADRD1$ox, q, na.rm = T)))
})] - (ox_plot[[2]]$fit[1])))



# GPS Matching
ADRD_pm25_match <- generate_pseudo_pop(ADRD1$N, ADRD1$pm25, 
                                       as.data.frame(ADRD1[, .(pct_male, pct_dual, pct_blk, pct_his, pct_oth_race, 
                                                 age, bmi, smoke, income, value, poverty,
                                                 ed, dens, ownocc, latitude, longitude)]),
                                       ci_appr = "matching",
                                       pred_model = "sl",
                                       gps_model = "parametric",
                                       use_cov_transform = TRUE,
                                       transformers = list("pow2", "pow3"),
                                       sl_lib = c("m_xgboost","SL.earth","SL.gam",
                                                  "SL.ranger"),
                                       params = list(xgb_nrounds=c(10, 20, 30),
                                                     xgb_eta=c(0.1, 0.2, 0.3)),
                                       nthread = 16,
                                       covar_bl_method = "absolute",
                                       covar_bl_trs = 0.1,
                                       trim_quantiles = c(0.01,0.99),
                                       max_attempt = 3,
                                       matching_fun = "matching_l1",
                                       delta_n = 1,
                                       scale = 0.5)
