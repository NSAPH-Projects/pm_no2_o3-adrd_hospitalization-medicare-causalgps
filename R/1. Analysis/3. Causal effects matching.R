###############################################################################
# Project: Causal effect of air pollution on first ADRD hospitalization       #
# Code:  #
# Input:                    
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
setDTthreads(threads = 8)

dir_data <- "/nfs/home/D/dam9096/shared_space/ci3_dam9096/Analysis_Data/Causal_ADRD/"

ADRDdat <- read_fst(paste0(dir_data, "ADRDhosp_all_counts_w_exposures.fst"), as.data.table = TRUE)
ADRD1 <- ADRDdat[, .(N = sum(n.cases),
                     size = sum(denom.size),
                     pct_male = sum((sex == 1) * denom.size) / sum(denom.size),
                     pct_dual = sum((dual == 1) * denom.size) / sum(denom.size),
                     pct_oth_race = sum((race %in% c(3, 4, 6)) * denom.size) / sum(denom.size),
                     pct_blk = sum((race == 2) * denom.size) / sum(denom.size),
                     pct_his = sum((race == 5) * denom.size) / sum(denom.size),
                     age = sum(avg_age * denom.size) / sum(denom.size),
                     bmi = mean(bmi),
                     smoke = mean(smoke),
                     # hispanic = mean(hispanic),
                     # pct_blk = mean(pct_blk),
                     income = mean(income),
                     value = mean(value),
                     poverty = mean(poverty),
                     ed = mean(education),
                     dens = mean(popdensity),
                     ownocc = mean(ownerocc),
                     pm25 = mean(pm25),
                     no2 = mean(no2),
                     ozone = mean(ozone),
                     ox = mean(ox)), 
                 by = c("year", "zip2")]
NE <- c("NY", "MA", "PA", "RI", "NH", "ME", "VT", "CT", "NJ")  
S <- c("DC", "VA", "NC", "WV", "KY", "SC", "GA", "FL", "AL", "TN", "MS", 
       "AR", "MD", "DE", "OK", "TX", "LA")
MW <- c("OH", "IN", "MI", "IA", "MO", "WI", "MN", "SD", "ND", "IL", "KS", "NE")
W <- c("MT", "CO", "WY", "ID", "UT", "NV", "CA", "OR", "WA", "AZ", "NM")
data(zipcode)
ADRD1 <- merge(ADRD1, zipcode, by.x = "zip2", by.y = "zip")
ADRD1[, region := ifelse(state %in% NE, "NE", ifelse(state %in% S, "S", ifelse(state %in% MW, "MW", "W")))]
ADRD1 <- ADRD1[complete.cases(ADRD1)]





# Estimate GPS
names(ADRD1)
ADRD_pm25_gps <- estimate_gps(ADRD1$N, ADRD1$pm25, 
                              ADRD1[, .(pct_male, pct_dual, pct_blk, pct_his, pct_oth_race, 
                                        age, bmi, smoke, income, value, poverty,
                                        ed, dens, ownocc, latitude, longitude)],
                              pred_model = "sl",
                              gps_model = "parametric",
                              internal_use = TRUE,
                              params = list(xgb_max_depth = c(3, 4, 5),
                                            xgb_rounds = c(10, 20, 30, 40)),
                              nthread = 16,                                
                              sl_lib = c("m_xgboost"))

# GLM
m <- glm(N ~ offset(log(size)) +
            pm25 + ox +
            region + as.factor(year) +
            poly(pct_dual, 2) + 
            pct_male + pct_oth_race + pct_blk + pct_his +
            poly(age, 2) + poly(bmi, 2) + 
            smoke + income + value +
            poverty + ed + dens + ownocc, 
          data = ADRD1, weights = ADRD_pm25_gps[[1]]$gps,
          family = poisson(link = "log"))
summary(m)

# GAM
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
