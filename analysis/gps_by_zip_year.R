rm(list = ls())
gc()

##### 0. Setup #####
# devtools::install_github("fasrc/CausalGPS", ref="develop")
library(data.table)
library(fst)
library(CausalGPS)
library(mgcv)
library(ggplot2)
library(tidyr)

# directory
dir_proj <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/"

# read in full data; 34,763,397 rows
ADRD_agg_lagged <- read_fst(paste0(dir_proj, "data/analysis/ADRD_complete_tv.fst"), as.data.table = TRUE)
setnames(ADRD_agg_lagged, old = c("pct_blk", "pct_owner_occ"), new = c("prop_blk", "prop_owner_occ"))
ADRD_agg_lagged[, `:=`(zip = as.factor(zip), year = as.factor(year), cohort = as.factor(cohort), age_grp = as.factor(age_grp), sex = as.factor(sex), race = as.factor(race), dual = as.factor(dual))]

source(paste0(dir_proj, "code/analysis/helper_functions.R"))

# parameters for this computing job
n_cores <- 16 # 48 or 64 if bigmem
n_gb <- 128 # 184 or 499 if bigmem
# total_n_rows <- nrow(ADRD_agg_lagged)
modifications <- "gps_by_zip_year" # to be used in names of output files, to record how you're tuning the models


##### Get data for exposure, outcome, and covariates of interest #####

exposure_name <- "pm25"
other_expos_names <- zip_expos_names[zip_expos_names != exposure_name]
outcome_name <- "n_hosp"
ADRD_agg_lagged_subset <- subset(ADRD_agg_lagged, select = c(exposure_name, outcome_name, other_expos_names, zip_var_names, "zip", indiv_var_names, offset_var_names))
for (var in c(zip_unordered_cat_var_names, indiv_unordered_cat_var_names)){
  ADRD_agg_lagged_subset[[var]] <- as.factor(ADRD_agg_lagged_subset[[var]])
}
colnames(ADRD_agg_lagged_subset)[colnames(ADRD_agg_lagged_subset) == exposure_name] <- "w"
# for (var in zip_ordered_cat_var_names){
#   ADRD_agg_lagged_subset[[var]] <- factor(ADRD_agg_lagged_subset[[var]], ordered = TRUE)
# }


##### Trim exposures outside the 1st and 99th percentiles, for all analyses (associational and causal) #####

trim_1_99 <- quantile(ADRD_agg_lagged_subset$w, c(0.01, 0.99))
rows_within_range <- ADRD_agg_lagged_subset$w >= trim_1_99[1] & ADRD_agg_lagged_subset$w <= trim_1_99[2]
ADRD_agg_lagged_trimmed_1_99 <- ADRD_agg_lagged_subset[rows_within_range, ]
n_rows <- nrow(ADRD_agg_lagged_trimmed_1_99) # 34,068,250


##### Estimate GPS for matching and weighting approaches #####

# get exposure and covariates, which are at ZIP, year level
zip_year_data <- subset(ADRD_agg_lagged_trimmed_1_99,
                        select = c("zip", "year", "w", other_expos_names, zip_var_names))
zip_year_data <- unique(zip_year_data, by = c("zip", "year"))
n_zip_year_rows <- nrow(zip_year_data) # 475,096

# create log file to see internal processes of CausalGPS
set_logger(logger_file_path = paste0(dir_proj, "code/analysis/CausalGPS_logs/CausalGPS_", Sys.Date(), "_estimate_gps_", modifications, "_", n_zip_year_rows, "rows_", n_cores, "cores_", n_gb, "gb.log"),
           logger_level = "TRACE")

# estimate GPS
set.seed(200)
zip_year_with_gps_obj <- estimate_gps(Y = 0, # fake Y variable since our outcomes are not at the zip-year level; not used in estimate_gps
                                  w = zip_year_data$w,
                                  c = subset(zip_year_data, select = c("year", other_expos_names, zip_var_names)),
                                  gps_model = "parametric", # i.e., model GPS normally, with mean and SD estimated using xgboost
                                  internal_use = T,
                                  params = list(xgb_nrounds = seq(10, 50),
                                                xgb_eta = seq(0.1, 0.4, 0.01)),
                                  sl_lib = c("m_xgboost"),
                                  nthread = n_cores)
zip_year_with_gps <- zip_year_with_gps_obj$dataset
zip_year_with_gps$zip <- zip_year_data$zip


##### GPS Weighting #####

# for GPS weighting, stabilize GPS using marginal probability of exposure (modeled normally) and cap extreme weights at 10
marginal_expos_prob <- dnorm(zip_year_data$w,
                             mean = mean(zip_year_data$w),
                             sd = sd(zip_year_data$w))
zip_year_with_gps$stabilized_ipw <- marginal_expos_prob / zip_year_with_gps$gps ## check estimate_gps
zip_year_with_gps$capped_stabilized_ipw <- ifelse(zip_year_with_gps$stabilized_ipw > 10, 10, zip_year_with_gps$stabilized_ipw)

# set up dataframe to check covariate balance
vars_for_cov_bal <- c("year", other_expos_names, zip_var_names)
abs_cor <- data.table(Covariate = rep(vars_for_cov_bal, 2),
                     Dataset = c(rep("Weighted", length(vars_for_cov_bal)), rep("Unweighted", length(vars_for_cov_bal))),
                     Absolute_Correlation = 100) # Absolute_Correlation column will be updated

# covariate balance: mean absolute point-biserial correlation for unordered categorical vars
for (unordered_var in c("year", zip_unordered_cat_var_names)){
  abs_cor[Covariate == unordered_var & Dataset == "Weighted", Absolute_Correlation := weighted_cor_unordered_var(zip_year_with_gps$w, zip_year_with_gps[[unordered_var]], zip_year_with_gps$capped_stabilized_ipw)]
  abs_cor[Covariate == unordered_var & Dataset == "Unweighted", Absolute_Correlation := cor_unordered_var(zip_year_with_gps$w, zip_year_with_gps[[unordered_var]])]
}

# covariate balance: absolute correlation for quantitative covariates
for (quant_var in c(other_expos_names, zip_quant_var_names)){
  abs_cor[Covariate == quant_var & Dataset == "Weighted", Absolute_Correlation := weightedCorr(zip_year_with_gps$w, zip_year_with_gps[[quant_var]], method = "Pearson", weights = zip_year_with_gps$capped_stabilized_ipw)]
  abs_cor[Covariate == quant_var & Dataset == "Unweighted", Absolute_Correlation := cor(zip_year_with_gps$w, zip_year_with_gps[[quant_var]])]
  abs_cor[Covariate == quant_var, Absolute_Correlation := abs(Absolute_Correlation)]
}

# plot covariate balance
weighted_cov_bal_plot <- ggplot(abs_cor, aes(x = Covariate, y = Absolute_Correlation, color = Dataset, group = Dataset)) +
  geom_point() +
  geom_line() +
  ylab(paste("Absolute Correlation with", exposure_name)) +
  ggtitle(paste("Set of", format(n_zip_year_rows, scientific = F, big.mark = ','), "ZIP-years")) +
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5))
ggsave(paste0(dir_proj, "results/covariate_balance/weighted_pop_", n_zip_year_rows, "rows", modifications, ".png"), weighted_cov_bal_plot)

# print summary statistics for pseudopopulation weights
ess(zip_year_with_gps$capped_stabilized_ipw) # 102,330.3

# if covariates balanced, model outcome from GPS-weighted pseudo-population
weighted_pseudopop <- merge(ADRD_agg_lagged_trimmed_1_99, subset(zip_year_with_gps,
                                                                 select = c("zip", "year", "capped_stabilized_ipw")),
                            by = c("zip", "year"))
indiv_var_names <- c("sex", "race", "dual", "age_grp")
formula_expos_only <- as.formula(paste(outcome_name, "~", paste(c("w", indiv_var_names), collapse = "+", sep = "")))
formula_expos_only_smooth <- as.formula(paste(outcome_name, "~", paste(c("s(w, bs = 'ts')", indiv_var_names), collapse = "+", sep = "")))

# parametric model (Poisson regression)
bam_exposure_only_capped_weighted <- bam(formula_expos_only,
                                         data = weighted_pseudopop,
                                         offset = log(n_persons * n_years),
                                         family = poisson(link = "log"),
                                         weights = capped_stabilized_ipw,
                                         samfrac = 0.05,
                                         chunk.size = 5000,
                                         control = gam.control(trace = TRUE))
summary(bam_exposure_only_capped_weighted)
saveRDS(summary(bam_exposure_only_capped_weighted), file = paste0(dir_proj, "results/parametric_results/bam_capped_weighted_exposure_only_", n_rows, "rows_", modifications, ".rds"))

# semi-parametric model (thin-plate spline)
bam_smooth_exposure_only_capped_weighted <- bam(formula_expos_only_smooth,
                                                data = weighted_pseudopop,
                                                offset = log(n_persons * n_years),
                                                family = poisson(link = "log"),
                                                weights = capped_stabilized_ipw,
                                                samfrac = 0.05,
                                                chunk.size = 5000,
                                                control = gam.control(trace = TRUE),
                                                nthreads = n_cores - 1)
png(paste0(dir_proj, "results/semiparametric_results/ERFs/bam_smooth_exposure_only_capped_weighted_", n_rows, "rows_", modifications, ".png"))
plot(bam_smooth_exposure_only_capped_weighted, main = paste0("GPS-Weighted, Capped at 10, Smoothed Poisson regression,\nexposure only (", exposure_name, ")"))
dev.off()
saveRDS(bam_smooth_exposure_only_capped_weighted, file = paste0(dir_proj, "results/semiparametric_results/spline_objects/bam_smooth_exposure_only_capped_weighted_", n_rows, "rows_", modifications, ".rds"))


##### GPS Matching #####
##### THE CODE BELOW IS UNFINISHED #####

# apply estimated GPS value to all strata within each ZIP/year
ADRD_agg_lagged_trimmed_1_99[, zip_year := .GRP, by = .(zip, year)]
zip_year_sizes <- as.vector(table(ADRD_agg_lagged_trimmed_1_99$zip_year))

expand_gps_values_to_units_of_analysis <- function(vector, gps_sizes){
  return(rep(vector, gps_sizes))
}
zip_year_with_gps_obj$e_gps_pred <- expand_gps_values_to_units_of_analysis(zip_year_with_gps_obj$e_gps_pred, zip_year_sizes)
zip_year_with_gps_obj$w_resid <- expand_gps_values_to_units_of_analysis(zip_year_with_gps_obj$w_resid, zip_year_sizes)

zip_year_with_gps_obj$dataset <- subset(zip_year_with_gps_obj$dataset,
                                        select = c("gps", "zip", "year"))
zip_year_with_gps_obj$dataset <- merge(zip_year_with_gps_obj$dataset, ADRD_agg_lagged_trimmed_1_99,
                                       by = c("zip", "year"))
zip_year_with_gps_obj$dataset[, `:=`(counter_weight = 0,
                                     row_index = 1:nrows,
                                     zip = NULL,
                                     zip_year = NULL,
                                     n_persons = NULL,
                                     n_years = NULL)]
# remove: indiv_var_names, "cohort", "year" ### check indiv_var_names (defined twice in this script)
# no actually, split into strata first? do all this within strata?
zip_year_with_gps_obj$dataset[[outcome_name]] <- NULL


# n_strata <- uniqueN(ADRD_agg_lagged_trimmed_1_99$stratum) # 19,714
# strata_sizes <- as.vector(table(ADRD_agg_lagged_trimmed_1_99$stratum))

##### THE CODE BELOW IS OLD #####

# GPS matching by ZIP-level covariates
set.seed(200)
matched_pop_subset <- generate_pseudo_pop(Y,
                                 w,
                                 c,
                                  ci_appr = "matching",
                                  pred_model = "sl",
                                  gps_model = "parametric",
                                  use_cov_transform = TRUE,
                                  transformers = list("sqrt", "log_nonneg", "logit_nonneg", "pow2", "pow3"),
                                  sl_lib = c("m_xgboost"),
                                  params = list(xgb_nrounds = seq(10, 50),
                                                xgb_eta = seq(0.1, 0.4, 0.01)),
                                  nthread = n_cores - 1,
                                  covar_bl_method = "absolute",
                                  covar_bl_trs = 0.1,
                                  covar_bl_trs_type = "maximal",
                                 optimized_compile = TRUE,
                                  trim_quantiles = c(0,1),
                                  max_attempt = 15,
                                  matching_fun = "matching_l1",
                                  delta_n = delta_n,
                                  scale = 1)
saveRDS(matched_pop_subset, file = paste0(dir_proj, "data/pseudopops/matched_pop_", n_rows, "rows", modifications, ".rds"))

# check ZIP-level covariate balance
# i.e., absolute correlation for quantitative covariates, polyserial correlation for ordered categorical variables, mean absolute point-biserial correlation for unordered categorical vars
matched_cov_bal_plot <- all_cov_bal(matched_pop_subset, w, c_unordered_vars = subset(c, select = zip_unordered_cat_var_names),
            ci_appr = "matching", all_cov_names = colnames(c), title = paste("Set of", format(n_rows, scientific = F), "observations"))
ggsave(paste0(dir_proj, "results/covariate_balance/matched_pop_", n_rows, "rows", modifications, ".png"), matched_cov_bal_plot)

# print summary statistics for pseudopopulation counter
summarize_pseudo_counter(matched_pop_subset)

# pseudopopulation, including individual-level covariates (i.e., strata), trimming away unmatched data
matched_obs <- matched_pop_subset$pseudo_pop$row_index[matched_pop_subset$pseudo_pop$counter > 0]
matched_indiv_vars <- indiv_vars[matched_obs] # To Do: consider including ffs_entry_year/ADRD_year in GPS or outcome model
matched_offset <- offset[matched_obs]
matched_data <- cbind(matched_pop_subset$pseudo_pop[matched_pop_subset$pseudo_pop$counter > 0, ], matched_indiv_vars, matched_offset) # to do: check that rows are in same order
matched_data <- as.data.frame(matched_data)

# Examine distribution of exposure and ZIP-level covariates (which were used to match) in pseudopopulation
summary(matched_data$w)
explore_zip_covs(matched_data)



##### IPTW by GPS using CausalGPS package #####

# create log file to see internal processes of CausalGPS
set_logger(logger_file_path = paste0(dir_proj, "code/analysis/CausalGPS_logs/CausalGPS_", Sys.Date(), "_weight_", modifications, "_", n_rows, "rows_", n_cores, "cores_", n_gb, "gb.log"),
           logger_level = "TRACE")

# GPS weighting on ZIP-level covariates
set.seed(200)
weighted_pop_subset <- generate_pseudo_pop(Y,
                                          w,
                                          c,
                                          ci_appr = "weighting",
                                          pred_model = "sl",
                                          gps_model = "parametric",
                                          use_cov_transform = TRUE,
                                          transformers = list("sqrt", "log_nonneg", "logit_nonneg", "pow2", "pow3"), # list("pow2", "pow3")
                                          sl_lib = c("m_xgboost"),
                                          params = list(xgb_nrounds = seq(10, 50),
                                                        xgb_eta = seq(0.1, 0.4, 0.01)),
                                          nthread = n_cores - 1,
                                          covar_bl_method = "absolute",
                                          covar_bl_trs = 0.2,
                                          covar_bl_trs_type = "maximal",
                                          optimized_compile = TRUE,
                                          trim_quantiles = c(0,1), # c(0.05, 0.95) or c(0.01, 0.99)
                                          max_attempt = 10,
                                          matching_fun = "matching_l1",
                                          delta_n = delta_n,
                                          scale = 1)
saveRDS(weighted_pop_subset, file = paste0(dir_proj, "data/pseudopops/weighted_pop_", n_rows, "rows", modifications, ".rds"))

# check ZIP-level covariate balance
# i.e., absolute correlation for quantitative covariates, polyserial correlation for ordered categorical variables, mean absolute point-biserial correlation for unordered categorical vars
weighted_cov_bal_plot <- all_cov_bal(weighted_pop_subset, w, c_unordered_vars = subset(c, select = zip_unordered_cat_var_names),
            "weighting", all_cov_names = colnames(c), title = paste("Set of", format(n_rows, scientific = F), "observations"))
ggsave(paste0(dir_proj, "results/covariate_balance/weighted_pop_", n_rows, "rows", modifications, ".png"), weighted_cov_bal_plot)

# print summary statistics for pseudopopulation weights
summarize_pseudo_weights(weighted_pop_subset)

# pseudopopulation, including individual-level covariates (i.e., strata), trimming away unmatched data
weighted_obs <- weighted_pop_subset$pseudo_pop$row_index
weighted_indiv_vars <- indiv_vars[weighted_obs] # To Do: consider including ffs_entry_year/ADRD_year in GPS or outcome model
weighted_offset <- offset[weighted_obs]
weighted_data <- cbind(weighted_pop_subset$pseudo_pop, weighted_indiv_vars, weighted_offset) # to do: check that rows are in same order
weighted_data <- as.data.frame(weighted_data)

# Examine distribution of exposure and ZIP-level covariates (which were used to match) in pseudopopulation
summary(weighted_data$w)
explore_zip_covs(weighted_data)


##### If desired, cap weights at 10 (or 99th percentile) #####

cap_weights = T
# cutoff_weight <- quantile(weighted_pop_subset$pseudo_pop$ipw, 0.99)
# print(cutoff_weight)
cutoff_weight <- 10

# to do: write function for this, taking cutoff as input
if (cap_weights){
  ipw <- weighted_pop_subset$pseudo_pop$counter_weight
  capped_weighted_pop_subset <- copy(weighted_pop_subset)
  capped_weighted_pop_subset$pseudo_pop$counter_weight <- ifelse(ipw > cutoff_weight, cutoff_weight, ipw)
  adjusted_corr_obj <- check_covar_balance(w = as.data.table(capped_weighted_pop_subset$pseudo_pop$w),
                                           c = subset(capped_weighted_pop_subset$pseudo_pop, select = c(other_expos_names, zip_quant_var_names)),
                                           ci_appr = "weighting",
                                           counter_weight = as.data.table(capped_weighted_pop_subset$pseudo_pop$counter_weight),
                                           nthread = n_cores - 1,
                                           covar_bl_method = "absolute",
                                           covar_bl_trs = 0.1, # or 0.1
                                           covar_bl_trs_type = "maximal",
                                           optimized_compile = T)
  capped_weighted_pop_subset$adjusted_corr_results <- adjusted_corr_obj$corr_results
}

cat("ESS of capped weighted pseudopopulation:", ess(capped_weighted_pop_subset$pseudo_pop$counter_weight))

# check ZIP-level covariate balance
# i.e., absolute correlation for quantitative covariates, polyserial correlation for ordered categorical variables, mean absolute point-biserial correlation for unordered categorical vars
capped_weighted_cov_bal_plot <- all_cov_bal(capped_weighted_pop_subset, w, c_unordered_vars = subset(c, select = zip_unordered_cat_var_names),
            "weighting", all_cov_names = colnames(c), title = paste("Set of", format(n_rows, scientific = F), "observations, weights capped at", cutoff_weight))
ggsave(paste0(dir_proj, "results/covariate_balance/capped_weighted_pop_", n_rows, "rows", modifications, ".png"), capped_weighted_cov_bal_plot)

# pseudopopulation, including individual-level covariates (i.e., strata), trimming away unmatched data
capped_weighted_obs <- capped_weighted_pop_subset$pseudo_pop$row_index
capped_weighted_indiv_vars <- indiv_vars[capped_weighted_obs] # To Do: consider including ffs_entry_year/ADRD_year in GPS or outcome model
capped_weighted_offset <- offset[capped_weighted_obs]
capped_weighted_data <- cbind(capped_weighted_pop_subset$pseudo_pop, capped_weighted_indiv_vars, capped_weighted_offset) # to do: check that rows are in same order
capped_weighted_data <- as.data.frame(capped_weighted_data)

# Examine distribution of exposure and ZIP-level covariates (which were used to match) in pseudopopulation
summary(capped_weighted_data$w)
explore_zip_covs(capped_weighted_data)


##### Old code: Estimate GPS using lm() #####

# estimate GPS using lm()
formula_lm_gps <- as.formula(paste("w ~", paste(c(other_expos_names, zip_var_names), collapse = "+", sep = "")))
lm_gps <- bam(formula_lm_gps, data = all_data)

# calculate raw GPS and IPW for data
lm_gps_data <- copy(all_data)
lm_gps_data[, gps := dnorm(w, lm_gps$fitted.values, sqrt(mean(lm_gps$residuals^2)))]
lm_gps_data[, ipw := 1/gps]

# stabilize GPS and IPW
w_mean <- mean(all_data$w)
w_sd <- sd(all_data$w)
lm_gps_data[, stab_gps := gps / dnorm(w, w_mean, w_sd)]
lm_gps_data[, stab_ipw := 1/stab_gps]
quantile(lm_gps_data$stab_ipw, c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1))

# truncate IPW at 10
lm_gps_data[, trunc_ipw := ifelse(stab_ipw > 10, 10, stab_ipw)]

# for UNTRUNCATED IPW, check absolute correlation for quantitative covariates # to do: other covariates
cor_val_pseudo <- sapply(subset(lm_gps_data, select = c(other_expos_names, zip_quant_var_names)), weightedCorr, y = lm_gps_data$w, method = "Pearson", weights = lm_gps_data$stab_ipw)
cor_val_orig <- sapply(subset(lm_gps_data, select = c(other_expos_names, zip_quant_var_names)), cor, lm_gps_data$w)
abs_cor = data.frame(Covariate = c(other_expos_names, zip_quant_var_names),
                     Unweighted = cor_val_orig,
                     Weighted = cor_val_pseudo) %>%
  gather(c(Unweighted, Weighted), key = 'Dataset', value = 'Absolute Correlation')
weighted_cov_bal_plot <- ggplot(abs_cor, aes(x = Covariate, y = `Absolute Correlation`, color = Dataset, group = Dataset)) +
  geom_point() +
  geom_line() +
  ggtitle(paste("Set of", format(n_rows, scientific = F), "observations, weights UNTRUNCATED")) +
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5))
ggsave(paste0(dir_proj, "results/covariate_balance/weighted_pop_", n_rows, "rows", modifications, ".png"), weighted_cov_bal_plot)

# for truncated IPW, check absolute correlation for quantitative covariates # to do: other covariates
cor_val_pseudo <- sapply(subset(lm_gps_data, select = c(other_expos_names, zip_quant_var_names)), weightedCorr, y = lm_gps_data$w, method = "Pearson", weights = lm_gps_data$trunc_ipw)
cor_val_orig <- sapply(subset(lm_gps_data, select = c(other_expos_names, zip_quant_var_names)), cor, lm_gps_data$w)
abs_cor = data.frame(Covariate = c(other_expos_names, zip_quant_var_names),
                     Unweighted = cor_val_orig,
                     Weighted = cor_val_pseudo) %>%
  gather(c(Unweighted, Weighted), key = 'Dataset', value = 'Absolute Correlation')
capped_weighted_cov_bal_plot <- ggplot(abs_cor, aes(x = Covariate, y = `Absolute Correlation`, color = Dataset, group = Dataset)) +
  geom_point() +
  geom_line() +
  ggtitle(paste("Set of", format(n_rows, scientific = F), "observations, weights truncated at 10")) +
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5))
ggsave(paste0(dir_proj, "results/covariate_balance/capped_weighted_pop_", n_rows, "rows", modifications, ".png"), capped_weighted_cov_bal_plot)


##### Old code #####

library(mgcv)

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