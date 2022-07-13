rm(list = ls())
gc()

##### 0. Setup #####
library(data.table)
library(fst)
# library(purrr)
# library(NSAPHutils)
library(CausalGPS)
library(mgcv)
library(ggplot2)
library(tidyr)

setDTthreads(threads = 16)
set.seed(100)

dir_proj <- "~/nsaph_projects/pm_no2_o3-adrd_hosp-medicare-causalgps/"
dir_data <- paste0(dir_proj, "data/")

ADRD_agg <- read_fst(paste0(dir_proj, "data/analysis/ADRD_complete_corrected.fst"), as.data.table = TRUE)

# Approximate first ADRD hospitalization by requiring no ADRD hosps for 2 years
ADRD_agg_lagged <- ADRD_agg[ADRD_year - ffs_entry_year >= 2, ]


#### Classify variables ####

exposure <- "pm25"
outcome <- "n_ADRDhosp"
zip_quant_var_names <- c("mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome",
                "medianhousevalue", "PIR", "poverty", "education", "popdensity", "pct_owner_occ",
                "summer_tmmx", "summer_rmax", "no2", "ozone_summer")
zip_var_names <- c(zip_quant_var_names, "ADRD_year", "region", "statecode")
zip_var_names2 <- c(zip_quant_var_names, "ADRD_year")
# zip_year_var_names <- c("ADRD_year", zip_var_names)

ADRD_agg_lagged_subset <- subset(ADRD_agg_lagged, select = c(exposure, outcome, zip_var_names))


#### Use CausalGPS package ####

# generate_pseudo_pop
# estimate_pmetric_erf
# estimate_npmetric_erf
# estimate_gps
# estimate_erf
# generate_syn_data


# dt_subset <- ADRD_agg_lagged[sexM == F & race_cat == "black" & any_dual == F & ADRD_age == 80 & ADRD_year == 2010 & ffs_entry_year == 2000]

# included_states <- "WY" # c("DC", "WY", "ND") # states with least observations
# dt_subset <- ADRD_agg_lagged_subset[statecode %in% included_states]
# dt_subset[, `:=`(region = NULL, statecode = NULL)]

# examine and/or transform quantitative covariates

# for (var in zip_quant_var_names){
#   hist(dt_subset[[var]], main = included_states, xlab = var)
# }
# skew_r_vars <- c("hispanic", "pct_blk", "medianhouseholdincome", "medianhousevalue", "PIR", "poverty", "popdensity")
# skew_l_vars <- c("pct_owner_occ")
# hist(sqrt(dt_subset$poverty+0.001))
# hist(sqrt(1-dt_subset$pct_owner_occ))

# Y_subset <- dt_subset$n_ADRDhosp
# w_subset <- dt_subset$pm25
# c_subset <- as.data.frame(subset(dt_subset, select = zip_var_names2))

# explore how many rows in c (i.e., exposures plus covariates) are duplicated
duplicates <- as.data.table(c_subset)[, .(count = .N), by = names(c_subset)]
summary(duplicates$count)

n_random_rows <- 35444
random_rows <- sample(1:nrow(ADRD_agg_lagged_subset), n_random_rows)
Y_subset <- ADRD_agg_lagged_subset[random_rows, n_ADRDhosp]
w_subset <- ADRD_agg_lagged_subset[random_rows, pm25]
c_subset <- as.data.frame(subset(ADRD_agg_lagged_subset[random_rows,], select = c("region", zip_var_names2)))

# Y <- ADRD_agg_lagged$n_ADRDhosp
# w <- ADRD_agg_lagged$pm25
# c <- as.data.frame(subset(ADRD_agg_lagged, select = zip_var_names)) # to do: consider including region, latitude, longitude


# To Do: transform covariates

# GPS matching on ZIP-level covariates
matched_pop_state <- generate_pseudo_pop(Y_subset,
                                 w_subset,
                                 c_subset,
                                  ci_appr = "matching",
                                  pred_model = "sl",
                                  gps_model = "parametric",
                                  use_cov_transform = TRUE,
                                  transformers = list("pow2", "pow3", "sqrt", "log"), # list("pow2", "pow3", "sqrt", "log")
                                  sl_lib = c("m_xgboost"),
                                  params = list(xgb_nrounds = c(10, 30, 50)),
                                  nthread = 15, # number of cores
                                  covar_bl_method = "absolute",
                                  covar_bl_trs = 0.15,
                                  covar_bl_trs_type = "maximal",
                                 optimized_compile = TRUE,
                                  trim_quantiles = c(0.05,0.95),
                                  max_attempt = 5,
                                  matching_fun = "matching_l1",
                                  delta_n = 0.1, # std dev of pm2.5 is 2.87, so I'll set delta_n = 0.2? parameters may depend on if state-level or national
                                  scale = 1)

# check ZIP-level covariate balance in matched data
cor_val_matched_state <- matched_pop_state$adjusted_corr_results
cor_val_unmatched_state <- matched_pop_state$original_corr_results
abs_cor = data.frame(cov = c(zip_quant_var_names, "region"), # colnames(c_subset)
                     unmatched = cor_val_unmatched_state$absolute_corr,
                     matched = cor_val_matched_state$absolute_corr) %>%
  gather(c(unmatched, matched), key = 'dataset', value = 'absolute correlation')
ggplot(abs_cor, aes(x = cov, y = `absolute correlation`, color = dataset, group = dataset)) +
  geom_point() +
  geom_line() +
  ggtitle(paste("Random set of", n_random_rows, "observations")) + # ggtitle(included_states)
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5))

# check number of matches
length(unique(matched_pop_state$pseudo_pop$row_index))
sum(matched_pop_state$pseudo_pop$counter)


##### IGNORE EVERYTHING BELOW THIS FOR NOW #####


# outcome model, including individual-level covariates (i.e., strata), trimming away unmatched data
# matched_obs <- matched_pop$pseudo_pop$row_index[matched_pop$pseudo_pop$counter > 0]
matched_obs <- matched_pop$pseudo_pop$row_index
indiv_vars <- ADRD_agg_lagged[matched_obs, .(sexM, race_cat, any_dual, ADRD_age)] # To Do: consider including ffs_entry_year/ADRD_year in GPS or outcome model
offset <- dt_subset[matched_obs, .(offset = n_persons * n_years)]
matched_data <- cbind(matched_pop$pseudo_pop, indiv_vars, offset) # to do: check that rows are in same order and whether this is data.frame

outcome <- estimate_pmetric_erf(formula = Y ~ . -row_index -gps -counter, # w or .?
                                  family = poisson, # poisson(link = "log")
                                  data = matched_data, # To Do: check if this must be a data.frame
                                  ci_appr = "matching")
summary(outcome)

# maybe just use gam instead of estimate_pmetric_erf

# non-linear model on matched data

erf_obj <- estimate_npmetric_erf(as.double(pseudo$Y),
                                 as.double(pseudo$w),
                                 bw_seq=seq(0.2,2,0.2),
                                 w_vals = seq(0,15,0.5),
                                 nthread = 16)

#summary(erf_obj)
plot(erf_obj)


#### Alternative: use gam ####



#### To Do: use transformed variables created in explore_covariate_distributions.R ####

# Estimate GPS using IPW

gps_pm25 <- bam(pm25 ~ factor(ADRD_year) + factor(any_dual) +
             ADRD_age + factor(sexM) + factor(race_cat) +
             mean_bmi + smoke_rate + hispanic + pct_blk +
             PIR + poverty + medhouseholdincome +
             education + popdensity + pct_owner_occ +
             no2 + ozone_summer,
           data = ADRD_agg_lagged,
           weights = ADRD_agg_lagged[, n_persons])

gps_no2 <- bam(no2 ~ factor(ADRD_year) + factor(any_dual) +
                 ADRD_age + factor(sexM) + factor(race_cat) +
                 mean_bmi + smoke_rate + hispanic + pct_blk +
                 PIR + poverty + medhouseholdincome +
                 education + popdensity + pct_owner_occ +
                 pm25 + ozone_summer,
               data = ADRD_agg_lagged,
               weights = ADRD_agg_lagged[, n_persons])

gps_ozone <- bam(ozone_summer ~ factor(ADRD_year) + factor(any_dual) +
                 ADRD_age + factor(sexM) + factor(race_cat) +
                 mean_bmi + smoke_rate + hispanic + pct_blk +
                 PIR + poverty + medhouseholdincome +
                 education + popdensity + pct_owner_occ +
                 pm25 + no2,
               data = ADRD_agg_lagged,
               weights = ADRD_agg_lagged[, n_persons])

# Stabilized IPWs
ADRD_agg_lagged[, `:=`(ipw_pm25 = dnorm(pm25, mean(pm25), sd(pm25)) /
           dnorm(pm25, gps_pm25$fitted.values, sqrt(mean(gps_pm25$residuals^2))))]
ADRD_agg_lagged[, `:=`(ipw_no2 = dnorm(no2, mean(no2), sd(no2)) /
             dnorm(no2, gps_no2$fitted.values, sqrt(mean(gps_no2$residuals^2))))]
ADRD_agg_lagged[, `:=`(ipw_ozone = dnorm(ozone_summer, mean(ozone_summer), sd(ozone_summer)) /
             dnorm(no2, gps_ozone$fitted.values, sqrt(mean(gps_ozone$residuals^2))))]
ADRD_agg_lagged[ipw_pm25 > 10, ipw_pm25 := 10]
ADRD_agg_lagged[ipw_no2 > 10, ipw_no2 := 10]
ADRD_agg_lagged[ipw_ozone > 10, ipw_ozone := 10]

# Covariate balance

bal <- cor(ADRD_agg_lagged[, pm25], 
           ADRD_agg_lagged[, .(ADRD_year, ADRD_age, sexM, any_dual,
                        white=race_cat=="white", black=race_cat=="black", other=race_cat=="other", 
                        hisp=race_cat=="hisp", n_am_nat=race_cat=="n_amer_native", asian=race_cat=="asian",
                        mean_bmi, smoke_rate, hispanic, pct_blk, PIR,
                        education, popdensity, pct_owner_occ, medhouseholdincome,
                        no2, ozone_summer)])
baldf <- data.frame(n = colnames(bal)[order(abs(bal))], 
                    y = 1:length(bal),
                    x = sort(abs(bal)))
baldf$n <- factor(baldf$n, levels = colnames(bal)[order(abs(bal))],
                  labels = colnames(bal)[order(abs(bal))])
bal_pm25 <- cov.wt(ADRD_agg_lagged[, .(pm25, 
                            ADRD_year, ADRD_age, sexM, any_dual,
                            white=race_cat=="white", black=race_cat=="black", 
                            other=race_cat=="other", hisp=race_cat=="hisp", 
                            n_am_nat=race_cat=="n_amer_native", asian=race_cat=="asian",
                            mean_bmi, smoke_rate, hispanic, pct_blk, PIR,
                            education, popdensity, pct_owner_occ, medhouseholdincome,
                            no2, ozone_summer)],
               wt = ADRD_agg_lagged[, ipw_pm25] * ADRD_agg_lagged[, n_persons],
               cor = TRUE)$cor[,1]
bal_no2 <- cov.wt(ADRD_agg_lagged[, .(no2, 
                                ADRD_year, ADRD_age, sexM, any_dual,
                                white=race_cat=="white", black=race_cat=="black", 
                                other=race_cat=="other", hisp=race_cat=="hisp", 
                                n_am_nat=race_cat=="n_amer_native", asian=race_cat=="asian",
                                mean_bmi, smoke_rate, hispanic, pct_blk, PIR,
                                education, popdensity, pct_owner_occ)],
                   wt = ADRD_agg_lagged[, ipw_no2] * ADRD_agg_lagged[, n_persons],
                   cor = TRUE)$cor[,1]
bal_ozone <- cov.wt(ADRD_agg_lagged[, .(ozone_summer, 
                                ADRD_year, ADRD_age, sexM, any_dual,
                                white=race_cat=="white", black=race_cat=="black", 
                                other=race_cat=="other", hisp=race_cat=="hisp", 
                                n_am_nat=race_cat=="n_amer_native", asian=race_cat=="asian",
                                mean_bmi, smoke_rate, hispanic, pct_blk, PIR,
                                education, popdensity, pct_owner_occ)],
                   wt = ADRD_agg_lagged[, ipw_ozone] * ADRD_agg_lagged[, n_persons],
                   cor = TRUE)$cor[,1]


bal2df <- merge(baldf, data.frame(n = names(bal_pm25), cor_pm25 = abs(bal_pm25)), by = c("n"))
bal2df <- merge(bal2df, data.frame(n = names(bal_no2), cor_no2 = abs(bal_no2)), by = c("n"))
bal2df <- merge(bal2df, data.frame(n = names(bal_ozone), cor_ozone = abs(bal_ozone)), by = c("n"))

ggplot(bal2df) +
  geom_point(aes(x = x, y = y)) +
  geom_line(aes(x = x, y = y)) +
  geom_point(aes(x = cor_pm25, y = y), color = 2) +
  # geom_point(aes(x = cor_no2, y = y), color = 3) +
  # geom_point(aes(x = cor_ozone, y = y), color = 4) +
  geom_vline(xintercept = 0.1, linetype = 2) +
  theme_bw() +
  scale_y_continuous(breaks = 1:length(bal), labels = baldf$n)


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

zip_weighted_pm25 <- unlist(map2(ADRD_agg_lagged$pm25, ADRD_agg_lagged$n_persons, rep))
zip_weighted_no2 <- unlist(map2(ADRD_agg_lagged$no2, ADRD_agg_lagged$n_persons, rep))
zip_weighted_ozone <- unlist(map2(ADRD_agg_lagged$ozone_summer, ADRD_agg_lagged$n_persons, rep))

cat("(Weighted by ZIP code aggregation) IQR of PM2.5 (micrograms/m^3):", IQR(ADRD_agg_lagged$pm25))
cat("(Weighted by ZIP code aggregation) IQR of NO2 (ppb):", IQR(ADRD_agg_lagged$no2))
cat("(Weighted by ZIP code aggregation) IQR of summer ozone (ppb):", IQR(ADRD_agg_lagged$ozone_summer))

cat("Hazard ratio per 1 IQR increase in PM2.5:", exp(model_pm25$coefficients[2]*IQR(zip_weighted_pm25)))
cat("Hazard ratio per 1 IQR increase in NO2:", exp(model_no2$coefficients[2]*IQR(zip_weighted_no2)))
cat("Hazard ratio per 1 IQR increase in summer ozone:", exp(model_ozone$coefficients[2]*IQR(zip_weighted_ozone)))

print("To compare with Shi et al 2021 preprint")
cat("Hazard ratio per 3.2 micrograms/m^3 increase in PM2.5:", exp(model_pm25$coefficients[2]*3.2))
cat("Hazard ratio per 11.6 ppb increase in NO2:", exp(model_no2$coefficients[2]*11.6))
cat("Hazard ratio per 5.3 ppb increase in summer ozone:", exp(model_ozone$coefficients[2]*5.3))

print("To compare with Shi et al 2020")
cat("Hazard ratio per 5 micrograms/m^3 increase in PM2.5:", exp(model_pm25$coefficients[2]*5))


##### Old code #####

# c <- copy(ADRD_agg_lagged)
# c[, `:=`(conf_year = NULL, year.y = NULL)] # conf_year = ADRD_year - 1, year.y is probably a duplicate variable
# c[, `:=`(n_ADRDhosp = NULL, pm25 = NULL, n_persons = NULL, n_years = NULL)]
# c[, `:=`(sexM = as.factor(sexM), race_cat = as.factor(race_cat), any_dual = as.factor(any_dual),
#          region = as.factor(region), # ADRD_zip = as.factor(ADRD_zip), city = as.factor(city), statecode = as.factor(statecode),
#          ADRD_year = as.factor(ADRD_year))] # ffs_entry_year = as.factor(ffs_entry_year)
# c[, `:=`(latitude = NULL, longitude = NULL, ADRD_zip = NULL, city = NULL, statecode = NULL,
#          ffs_entry_year = NULL)]

# model.matrix(y~-1+as.factor(year)+age,data)
# c <- as.matrix(c) # don't use this
# c <- as.data.frame(c)