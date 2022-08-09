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

# To Do: consider using more cores
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
offset_var_names <- c("n_persons", "n_years")

zip_quant_var_names <- c("mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome",
                "medianhousevalue", "PIR", "poverty", "education", "popdensity", "pct_owner_occ",
                "summer_tmmx", "summer_rmax", "no2", "ozone_summer")
zip_unordered_cat_var_names <- c("region", "ADRD_year")

indiv_quant_var_names <- c("ADRD_age")
indiv_unordered_cat_var_names <- c("sexM", "race_cat", "any_dual")

zip_var_names <- c(zip_quant_var_names, zip_unordered_cat_var_names)
indiv_var_names <- c(indiv_unordered_cat_var_names, indiv_quant_var_names) # note: for now, using ADRD_age as a quantitative variable (not binned)

ADRD_agg_lagged_subset <- subset(ADRD_agg_lagged, select = c(exposure_name, outcome_name, zip_var_names, indiv_var_names, offset_var_names))
for (var in c(zip_unordered_cat_var_names, indiv_unordered_cat_var_names)){
  ADRD_agg_lagged_subset[[var]] <- as.factor(ADRD_agg_lagged_subset[[var]])
}
# for (var in zip_ordered_cat_var_names){
#   ADRD_agg_lagged_subset[[var]] <- factor(ADRD_agg_lagged_subset[[var]], ordered = TRUE)
# }


##### Subset data to make code run faster #####

set.seed(100)
n_random_rows <- 100000 # To Do: use a bigger sample
random_rows <- sample(1:nrow(ADRD_agg_lagged_subset), n_random_rows)
Y_subset <- ADRD_agg_lagged_subset[random_rows, n_ADRDhosp]
w_subset <- ADRD_agg_lagged_subset[random_rows, pm25]
c_subset <- as.data.frame(subset(ADRD_agg_lagged_subset[random_rows,], select = zip_var_names))

# To Do at later stage: Full data
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

# To Do: consider using larger delta_n, more cores, lm [or another package Kevin may have mentioned?] instead of xgboost, more attempts

# GPS matching by ZIP-level covariates
set.seed(200)
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
                                  nthread = 15,
                                  covar_bl_method = "absolute",
                                  covar_bl_trs = 0.1,
                                  covar_bl_trs_type = "maximal",
                                 optimized_compile = TRUE,
                                  trim_quantiles = c(0.05,0.95),
                                  max_attempt = 10,
                                  matching_fun = "matching_l1",
                                  delta_n = 0.1,
                                  scale = 1)

# check ZIP-level covariate balance in matched data: abs correlation for quantitative or ordered categorical variables
cor_val_unmatched_subset <- matched_pop_subset$original_corr_results$absolute_corr[zip_quant_var_names] # remove non-ordinal categorical variables; can include zip_ordered_cat_var_names if exists
cor_val_matched_subset <- matched_pop_subset$adjusted_corr_results$absolute_corr[zip_quant_var_names]
abs_cor = data.frame(Covariate = zip_quant_var_names,
                     Unmatched = cor_val_unmatched_subset,
                     Matched = cor_val_matched_subset) %>%
  gather(c(Unmatched, Matched), key = 'Dataset', value = 'Absolute Correlation')
ggplot(abs_cor, aes(x = Covariate, y = `Absolute Correlation`, color = Dataset, group = Dataset)) +
  geom_point() +
  geom_line() +
  ggtitle(paste("Random set of", n_random_rows, "observations")) + # ggtitle(included_states)
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5))

# check ZIP-level covariate balance in matched data: unordered categorical variables
# to do: contrast test or something, consider writing function or for loop for all unordered cat vars
ggplot(matched_pop_subset$pseudo_pop, aes(x = region, y = w, weight = counter)) +
  geom_boxplot()
ggplot(matched_pop_subset$pseudo_pop, aes(x = ADRD_year, y = w, weight = counter)) +
  geom_boxplot()

# check number of matches
counter <- matched_pop_subset$pseudo_pop$counter
cat("Number of observations included in pseudo-population:", length(unique(matched_pop_subset$pseudo_pop$row_index)))
cat("Total number of matches:", sum(counter))
print("Distribution of number of matches per included observation:")
summary(counter)

# check Kish's effective sample size
ess <- sum(counter)^2 / (sum(counter^2))
ess


##### If desired, cap matches ("counter" variable) at 95th percentile #####
# July 27 edit: don't do this for matching

cap_counts = F

if (cap_counts){
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
  cor_val_unmatched_subset <- matched_pop_subset$original_corr_results$absolute_corr[zip_quant_var_names] # remove non-ordinal categorical variables; can include zip_ordered_cat_var_names if exists
  cor_val_matched_subset <- matched_pop_subset$adjusted_corr_results$absolute_corr[zip_quant_var_names]
  abs_cor = data.frame(Covariate = zip_quant_var_names,
                       Unmatched = cor_val_unmatched_subset,
                       Matched = cor_val_matched_subset) %>%
    gather(c(Unmatched, Matched), key = 'Dataset', value = 'Absolute Correlation')
  ggplot(abs_cor, aes(x = Covariate, y = `Absolute Correlation`, color = Dataset, group = Dataset)) +
    geom_point() +
    geom_line() +
    ggtitle(paste("Random set of", n_random_rows, "observations")) + # ggtitle(included_states)
    theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5))
}


##### IPTW by GPS using CausalGPS package #####

# To Do: consider using larger delta_n, more cores, lm [or another package Kevin may have mentioned?] instead of xgboost, more attempts

# GPS weighting on ZIP-level covariates
set.seed(200)
weighted_pop_subset <- generate_pseudo_pop(Y_subset,
                                          w_subset,
                                          c_subset,
                                          ci_appr = "weighting",
                                          pred_model = "sl",
                                          gps_model = "parametric",
                                          use_cov_transform = TRUE,
                                          transformers = list("pow2", "pow3", "sqrt", "log_nonneg", "logit_nonneg"), # list("pow2", "pow3")
                                          sl_lib = c("m_xgboost"),
                                          params = list(xgb_nrounds = c(10, 20, 30, 50)),
                                          nthread = 15,
                                          covar_bl_method = "absolute",
                                          covar_bl_trs = 0.1,
                                          covar_bl_trs_type = "maximal",
                                          optimized_compile = TRUE,
                                          trim_quantiles = c(0.05,0.95), # c(0.05, 0.95) or c(0.01, 0.99)
                                          max_attempt = 10,
                                          matching_fun = "matching_l1",
                                          delta_n = 0.1, # std dev of pm2.5 is 2.87, so I'll set delta_n = 0.2? parameters may depend on if state-level or national
                                          scale = 1)

# check ZIP-level covariate balance in matched data: abs correlation for quantitative or ordered categorical variables
cor_val_unweighted_subset <- weighted_pop_subset$original_corr_results$absolute_corr[zip_quant_var_names] # remove non-ordinal categorical variables; can include zip_ordered_cat_var_names if exists
cor_val_weighted_subset <- weighted_pop_subset$adjusted_corr_results$absolute_corr[zip_quant_var_names]
abs_cor_weighting = data.frame(Covariate = zip_quant_var_names,
                     Unweighted = cor_val_unweighted_subset,
                     Weighted = cor_val_weighted_subset) %>%
  gather(c(Unweighted, Weighted), key = 'Dataset', value = 'Absolute Correlation')
ggplot(abs_cor_weighting, aes(x = Covariate, y = `Absolute Correlation`, color = Dataset, group = Dataset)) +
  geom_point() +
  geom_line() +
  ggtitle(paste("Random set of", n_random_rows, "observations")) + # ggtitle(included_states)
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5))

# check ZIP-level covariate balance in matched data: unordered categorical variables
# to do: contrast test or something, consider writing function or for loop for all unordered cat vars
ggplot(weighted_pop_subset$pseudo_pop, aes(x = region, y = w, weight = ipw)) +
  geom_boxplot()
ggplot(weighted_pop_subset$pseudo_pop, aes(x = ADRD_year, y = w, weight = ipw)) +
  geom_boxplot()

# check number of matches
ipw <- weighted_pop_subset$pseudo_pop$ipw
cat("Number of observations included in pseudo-population:", length(unique(matched_pop_subset$pseudo_pop$row_index)))
cat("Sum of weights:", sum(ipw))
print("Distribution of weights:")
summary(ipw)
quantile(ipw, c(0, 0.25, 0.5, 0.75, 0.95, 0.99, 0.999))

# check Kish's effective sample size
ipw <- weighted_pop_subset$pseudo_pop$ipw
ess <- sum(ipw)^2 / (sum(ipw^2))
ess


##### If desired, cap weights at 95th percentile #####

cap_weights = T

if (cap_weights){
  weighted_pop_subset_wins <- copy(weighted_pop_subset)
  cutoff_weight <- quantile(weighted_pop_subset_wins$pseudo_pop$ipw, 0.95)
  weighted_pop_subset_wins$pseudo_pop$ipw <- ifelse(weighted_pop_subset_wins$pseudo_pop$ipw > cutoff_weight, cutoff_weight, weighted_pop_subset_wins$pseudo_pop$ipw)
  adjusted_corr_obj <- check_covar_balance(weighted_pop_subset_wins$pseudo_pop,
                                           ci_appr="weighting",
                                           nthread=15,
                                           covar_bl_method = "absolute",
                                           covar_bl_trs = 0.1,
                                           covar_bl_trs_type = "maximal",
                                           optimized_compile=T)
  weighted_pop_subset_wins$adjusted_corr_results <-  adjusted_corr_obj$corr_results
  
  # check ZIP-level covariate balance in matched data: abs correlation for quantitative or ordered categorical variables
  cor_val_unweighted_subset <- weighted_pop_subset_wins$original_corr_results$absolute_corr[zip_quant_var_names] # remove non-ordinal categorical variables; can include zip_ordered_cat_var_names if exists
  cor_val_weighted_subset <- weighted_pop_subset_wins$adjusted_corr_results$absolute_corr[zip_quant_var_names]
  abs_cor_weighting = data.frame(Covariate = zip_quant_var_names,
                                 Unweighted = cor_val_unweighted_subset,
                                 Weighted = cor_val_weighted_subset) %>%
    gather(c(Unweighted, Weighted), key = 'Dataset', value = 'Absolute Correlation')
  ggplot(abs_cor_weighting, aes(x = Covariate, y = `Absolute Correlation`, color = Dataset, group = Dataset)) +
    geom_point() +
    geom_line() +
    ggtitle(paste("Random set of", n_random_rows, "observations")) + # ggtitle(included_states)
    theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5))
}


##### (Poisson) Parametric outcome models #####


# outcome model, including individual-level covariates (i.e., strata), trimming away unmatched data
matched_obs <- matched_pop_subset$pseudo_pop$row_index[matched_pop_subset$pseudo_pop$counter > 0]
matched_indiv_vars <- indiv_vars_subset[matched_obs] # To Do: consider including ffs_entry_year/ADRD_year in GPS or outcome model
matched_offset <- offset_subset[matched_obs]
matched_data <- cbind(matched_pop_subset$pseudo_pop[matched_pop_subset$pseudo_pop$counter > 0, ], matched_indiv_vars, matched_offset) # to do: check that rows are in same order
matched_data <- as.data.frame(matched_data)

# weighted model
weighted_obs <- matched_pop_subset$pseudo_pop$row_index
weighted_indiv_vars <- indiv_vars_subset[weighted_obs] # To Do: consider including ffs_entry_year/ADRD_year in GPS or outcome model
weighted_offset <- offset_subset[weighted_obs]
weighted_data <- cbind(weighted_pop_subset$pseudo_pop, weighted_indiv_vars, weighted_offset) # to do: check that rows are in same order
weighted_data <- as.data.frame(weighted_data)

# method 1: gam package
# To do: rd rest of Kevin's code (https://github.com/kevjosey/erc-strata/blob/master/R/match_estimate.R)
# to see if there's anything extra I should do to account for matching

# GPS matching
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

bam_exposure_only_matched <- bam(Y ~ w + any_dual + ADRD_age + sexM + race_cat,
                         data = matched_data,
                         offset = log(person_years),
                         family = poisson(link = "log"),
                         weights = counter,
                         samfrac = 0.05,
                         chunk.size = 5000,
                         control = gam.control(trace = TRUE))
summary(bam_exposure_only_matched)

# GPS weighting
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

bam_exposure_only_weighted <- bam(Y ~ w + any_dual + ADRD_age + sexM + race_cat,
                                 data = weighted_data,
                                 offset = log(person_years),
                                 family = poisson(link = "log"),
                                 weights = ipw,
                                 samfrac = 0.05,
                                 chunk.size = 5000,
                                 control = gam.control(trace = TRUE))
summary(bam_exposure_only_weighted)


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
semipmetric_exposure_only <- estimate_semipmetric_erf(Y ~ w + any_dual + ADRD_age + sexM + race_cat,
                                                      data = matched_data,
                                                      # offset = log(person_years),
                                                      family = poisson(link = "log"),
                                                      ci_appr = "matching")
plot(semipmetric_exposure_only)



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


# model.matrix(y~-1+as.factor(year)+age,data)
# c <- as.matrix(c) # don't use this
# c <- as.data.frame(c)