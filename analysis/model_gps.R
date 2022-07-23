rm(list = ls())
gc()

##### 0. Setup #####
library(data.table)
library(fst)
# library(purrr)
# library(NSAPHutils)
library(CausalGPS)
library(mgcv)
library(gnm)
library(ggplot2)
library(tidyr)

setDTthreads(threads = 16)
set.seed(100)

dir_proj <- "~/nsaph_projects/pm_no2_o3-adrd_hosp-medicare-causalgps/"
dir_data <- paste0(dir_proj, "data/")

# create log file to see internal processes of CausalGPS
# to do: change name of repo from "pm_no2_o3-adrd_hospitalization-medicare-causalgps" (since name is similar to project folder) to "git"
set_logger(logger_file_path = paste0(dir_proj, "pm_no2_o3-adrd_hospitalization-medicare-causalgps/analysis/CausalGPS.log"),
           logger_level = "TRACE")

ADRD_agg <- read_fst(paste0(dir_proj, "data/analysis/ADRD_complete_corrected.fst"), as.data.table = TRUE)

# Approximate first ADRD hospitalization by requiring no ADRD hosps for 2 years
ADRD_agg_lagged <- ADRD_agg[ADRD_year - ffs_entry_year >= 2, ]


##### Classify variables #####

exposure_name <- "pm25"
outcome_name <- "n_ADRDhosp"
zip_quant_var_names <- c("mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome",
                "medianhousevalue", "PIR", "poverty", "education", "popdensity", "pct_owner_occ",
                "summer_tmmx", "summer_rmax", "no2", "ozone_summer")
zip_cat_var_names <- c("region", "ADRD_year")
zip_var_names <- c(zip_quant_var_names, zip_cat_var_names)
indiv_cat_var_names <- c("sexM", "race_cat", "any_dual")
indiv_var_names <- c(indiv_cat_var_names, "ADRD_age") # note: for now, using ADRD_age as a quantitative variable (not binned)
offset_var_names <- c("n_persons", "n_years")

ADRD_agg_lagged_subset <- subset(ADRD_agg_lagged, select = c(exposure_name, outcome_name, zip_var_names, indiv_var_names, offset_var_names))
for (var in c(zip_cat_var_names, indiv_cat_var_names)){
  ADRD_agg_lagged_subset[[var]] <- as.factor(ADRD_agg_lagged_subset[[var]])
}


##### Subset data to make code run faster #####

set.seed(100)
n_random_rows <- 100000
random_rows <- sample(1:nrow(ADRD_agg_lagged_subset), n_random_rows)
# Y_subset <- subset(ADRD_agg_lagged_subset, subset = (row.names(ADRD_agg_lagged_subset) == random_rows), select = outcome_name)
# w_subset <- subset(ADRD_agg_lagged_subset, subset = random_rows, select = exposure_name)
Y_subset <- ADRD_agg_lagged_subset[random_rows, n_ADRDhosp]
w_subset <- ADRD_agg_lagged_subset[random_rows, pm25]
c_subset <- as.data.frame(subset(ADRD_agg_lagged_subset[random_rows,], select = zip_var_names))

# Full data
# Y <- ADRD_agg_lagged$n_ADRDhosp
# w <- ADRD_agg_lagged$pm25
# c <- as.data.frame(subset(ADRD_agg_lagged, select = zip_var_names))

# Not used in GPS matching, but used in outcome model
indiv_vars_subset <- ADRD_agg_lagged_subset[random_rows, .(sexM, race_cat, any_dual, ADRD_age)] # To Do: consider including ffs_entry_year/ADRD_year in GPS or outcome model
offset_subset <- ADRD_agg_lagged_subset[random_rows, .(person_years = n_persons * n_years)]
data_subset <- cbind(Y_subset, w_subset, indiv_vars_subset, c_subset, offset_subset)


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


##### Write functions to transform variables #####

log_nonneg <- function(x) log(x + 0.001)
logit_nonneg <- function(x) log((x + 0.001)/(1 - x + 0.001))


##### Match on GPS using CausalGPS package #####

# GPS matching on ZIP-level covariates
matched_pop_subset <- generate_pseudo_pop(Y_subset,
                                 w_subset,
                                 c_subset,
                                  ci_appr = "matching",
                                  pred_model = "sl",
                                  gps_model = "parametric",
                                  use_cov_transform = TRUE,
                                  transformers = list("pow2", "pow3", "sqrt", "log_nonneg", "logit_nonneg"), # list("pow2", "pow3")
                                  sl_lib = c("m_xgboost"),
                                  params = list(xgb_nrounds = c(10, 20, 30, 50)),
                                  nthread = 15, # number of cores
                                  covar_bl_method = "absolute",
                                  covar_bl_trs = 0.1,
                                  covar_bl_trs_type = "maximal",
                                 optimized_compile = TRUE,
                                  trim_quantiles = c(0.05,0.95), # c(0.05, 0.95) or c(0.01, 0.99)
                                  max_attempt = 5,
                                  matching_fun = "matching_l1",
                                  delta_n = 0.1, # std dev of pm2.5 is 2.87, so I'll set delta_n = 0.2? parameters may depend on if state-level or national
                                  scale = 1)

# check ZIP-level covariate balance in matched data
cor_val_matched_subset <- matched_pop_subset$adjusted_corr_results
cor_val_unmatched_subset <- matched_pop_subset$original_corr_results
abs_cor = data.frame(covariate = zip_var_names, # colnames(c_subset)
                     unmatched = cor_val_unmatched_subset$absolute_corr,
                     matched = cor_val_matched_subset$absolute_corr) %>%
  gather(c(unmatched, matched), key = 'dataset', value = 'absolute correlation')
ggplot(abs_cor, aes(x = covariate, y = `absolute correlation`, color = dataset, group = dataset)) +
  geom_point() +
  geom_line() +
  ggtitle(paste("Random set of", n_random_rows, "observations")) + # ggtitle(included_states)
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5))

# check number of matches
counter <- matched_pop_subset$pseudo_pop$counter
cat("Number of observations included in pseudo-population:", length(unique(matched_pop_subset$pseudo_pop$row_index)))
cat("Total number of matches:", sum(counter))
print("Distribution of number of matches per included observation:")
summary(counter)

# To Do: fix this; ess1 and ess2 are not equivalent; both formulas could be wrong
# check effective sample size
# source of formula: https://stats.stackexchange.com/questions/499397/matchit-output-using-coarsened-exact-matching
ess1 <- sum(counter)^2 / (sum(counter^2))
counter_scaled_to_1 <- counter/sum(counter) # scale weights to sum to 1
ess2 <- length(counter_scaled_to_1) / (1 + var(counter_scaled_to_1)) # should be equivalent to ess1
ess1
ess2


##### If desired, winsorize matches ("counter" variable) at 95th percentile #####

winsorize_counts = F

if (winsorize_counts){
  cutoff <- quantile(matched_pop_subset$pseudo_pop$counter, 0.95)
  matched_pop_subset$pseudo_pop$counter <- ifelse(matched_pop_subset$pseudo_pop$counter > cutoff, cutoff, matched_pop_subset$pseudo_pop$counter)
  adjusted_corr_obj <- check_covar_balance(matched_pop_subset$pseudo_pop,
                                           ci_appr="matching",
                                           nthread=15,
                                           covar_bl_method = "absolute",
                                           covar_bl_trs = 0.1,
                                           covar_bl_trs_type = "maximal",
                                           optimized_compile=T)
  matched_pop_subset$adjusted_corr_results <-  adjusted_corr_obj$corr_results
  
  # check ZIP-level covariate balance in matched data
  cor_val_matched_subset <- matched_pop_subset$adjusted_corr_results
  cor_val_unmatched_subset <- matched_pop_subset$original_corr_results
  abs_cor = data.frame(covariate = zip_var_names, # colnames(c_subset)
                       unmatched = cor_val_unmatched_subset$absolute_corr,
                       matched = cor_val_matched_subset$absolute_corr) %>%
    gather(c(unmatched, matched), key = 'dataset', value = 'absolute correlation')
  ggplot(abs_cor, aes(x = covariate, y = `absolute correlation`, color = dataset, group = dataset)) +
    geom_point() +
    geom_line() +
    ggtitle(paste("Random set of", n_random_rows, "observations")) + # ggtitle(included_states)
    theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5))
}


##### (Poisson) Parametric outcome models #####


# outcome model, including individual-level covariates (i.e., strata), trimming away unmatched data
# matched_obs <- matched_pop_subset$pseudo_pop$row_index[matched_pop_subset$pseudo_pop$counter > 0]
matched_obs <- matched_pop_subset$pseudo_pop$row_index
matched_indiv_vars <- indiv_vars_subset[matched_obs] # To Do: consider including ffs_entry_year/ADRD_year in GPS or outcome model
matched_offset <- offset_subset[matched_obs]
matched_data <- cbind(matched_pop_subset$pseudo_pop, matched_indiv_vars, matched_offset) # to do: check that rows are in same order
matched_data <- as.data.frame(matched_data)

# method 1: gam package
# To do: figure out if exposure_only and exposures_controlled models should stratify across indiv characteristics
# To do: rd rest of Kevin's code to see if there's anything extra I should do to account for matching
bam_doubly_robust <- bam(Y ~ w + no2 + ozone_summer +
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
summary(bam_doubly_robust)

bam_exposures_controlled <- bam(Y ~ w + no2 + ozone_summer,
                                data = matched_data,
                                offset = log(person_years),
                                family = poisson(link = "log"),
                                weights = counter,
                                samfrac = 0.05,
                                chunk.size = 5000,
                                control = gam.control(trace = TRUE))
summary(bam_exposures_controlled)

bam_exposure_only <- bam(Y ~ w,
                    data = matched_data,
                    offset = log(person_years),
                    family = poisson(link = "log"),
                    weights = counter,
                    samfrac = 0.05,
                    chunk.size = 5000,
                    control = gam.control(trace = TRUE))
summary(bam_exposure_only)


# method 2: gnm package
# To do: figure out how to use eliminate (rn, there's an error: only first element used)
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
# To do: wrong right now cuz offset shouldn't be constant; to do - fix
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


##### Non-parametric model on matched data ##### 

erf_obj <- estimate_npmetric_erf(as.double(matched_data$Y),
                                 as.double(matched_data$w),
                                 matched_data$counter, # if I comment out this line, the ERF is much smoother (why?)
                                 bw_seq=seq(0.2,2,0.2),
                                 w_vals = seq(0,15,0.5),
                                 nthread = 16)
plot(erf_obj)
#summary(erf_obj)



#### To Do: use transformed variables created in explore_covariate_distributions.R ####


##### Old code #####

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