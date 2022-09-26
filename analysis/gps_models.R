rm(list = ls())
gc()

##### 0. Setup #####
library(data.table)
library(fst)
library(CausalGPS)
library(mgcv)
# library(gnm)
library(ggplot2)
library(tidyr)

# directory
dir_proj <- "~/nsaph_projects/pm_no2_o3-adrd_hosp-medicare-causalgps/"
# dir_data <- paste0(dir_proj, "data/")
# dir_code <- paste0(dir_proj, "code/")

# parameters for this computing job
n_cores <- 48
n_gb <- 184
# n_rows <- 5000000
n_rows <- 17640610 # Full data
modifications <- "" # to be used in names of output files, for instance "delta0.3_"


# read in full data
ADRD_agg <- read_fst(paste0(dir_proj, "data/analysis/ADRD_complete_corrected.fst"), as.data.table = TRUE)
ADRD_agg_lagged <- ADRD_agg[ADRD_year - ffs_entry_year >= 2, ] # Approximate first ADRD hospitalization by requiring no ADRD hosps for 2 years
# bin ages into 5-year increments (except first/last bin are smaller/larger)
# ADRD_agg_lagged[, ADRD_age := ifelse(ADRD_age < 70, 1, ifelse(ADRD_age < 75, 2, ifelse(ADRD_age < 80, 3, ifelse(ADRD_age < 85, 4, ifelse(ADRD_age < 90, 5, ifelse(ADRD_age < 95, 6, 7))))))]

source(paste0(dir_proj, "code/analysis/helper_functions.R"))


##### Get data for exposure, outcome, and covariates of interest #####

exposure_name <- "pm25"
outcome_name <- "n_ADRDhosp"
ADRD_agg_lagged_subset <- subset(ADRD_agg_lagged, select = c(exposure_name, outcome_name, zip_var_names, indiv_var_names, offset_var_names))
for (var in c(zip_unordered_cat_var_names, indiv_unordered_cat_var_names)){
  ADRD_agg_lagged_subset[[var]] <- as.factor(ADRD_agg_lagged_subset[[var]])
}
# for (var in zip_ordered_cat_var_names){
#   ADRD_agg_lagged_subset[[var]] <- factor(ADRD_agg_lagged_subset[[var]], ordered = TRUE)
# }


##### Use full data or (random) subset of data to make code run faster than full data #####

if (n_rows < 17640610){ # if analyzing subset of data
  set.seed(100)
  selected_rows <- sample(1:nrow(ADRD_agg_lagged_subset), n_rows)
} else selected_rows <- 1:17640610 # if full data

Y_subset <- ADRD_agg_lagged_subset[selected_rows, n_ADRDhosp]
w_subset <- ADRD_agg_lagged_subset[selected_rows, pm25]
c_subset <- as.data.frame(subset(ADRD_agg_lagged_subset[selected_rows,], select = zip_var_names))

# Not used in GPS matching, but used in outcome model
indiv_vars_subset <- ADRD_agg_lagged_subset[selected_rows, .(sexM, race_cat, any_dual, ADRD_age)] # To Do: consider including ffs_entry_year/ADRD_year in GPS or outcome model
offset_subset <- ADRD_agg_lagged_subset[selected_rows, .(person_years = n_persons * n_years)]
data_subset <- cbind(Y_subset, w_subset, indiv_vars_subset, c_subset, offset_subset)

# Examine distribution of exposure and ZIP-level covariates in sample
summary(w_subset)
explore_zip_covs(c_subset)

# To Do at later stage: Full data
# Y <- ADRD_agg_lagged$n_ADRDhosp
# w <- ADRD_agg_lagged$pm25
# c <- as.data.frame(subset(ADRD_agg_lagged, select = zip_var_names))


##### Match on GPS using CausalGPS package #####

# create log file to see internal processes of CausalGPS
set_logger(logger_file_path = paste0(dir_proj, "code/analysis/CausalGPS_logs/CausalGPS_", Sys.Date(), "_match_", modifications, n_rows, "rows_", n_cores, "cores_", n_gb, "gb.log"),
           logger_level = "DEBUG")

# if using SL.gam, remove mgcv library and allow custom parameters
# library(SuperLearner)
# m_gam <- function(cts.num = 16, deg.gam = 1, ...) SL.gam(cts.num = cts.num, deg.gam = deg.gam, ...) # to do: explain these params
# detach("package:mgcv", unload=TRUE)

# GPS matching by ZIP-level covariates
set.seed(200)
matched_pop_subset <- generate_pseudo_pop(Y_subset,
                                 w_subset,
                                 c_subset,
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
                                  covar_bl_trs = 0.2,
                                  covar_bl_trs_type = "maximal",
                                 optimized_compile = TRUE,
                                  trim_quantiles = c(0.01,0.99),
                                  max_attempt = 10,
                                  matching_fun = "matching_l1",
                                  delta_n = 0.2,
                                  scale = 1)
saveRDS(matched_pop_subset, file = paste0(dir_proj, "data/pseudopops/matched_pop_", n_rows, "rows", modifications, ".rds"))

# check ZIP-level covariate balance
# i.e., absolute correlation for quantitative covariates, polyserial correlation for ordered categorical variables, mean absolute point-biserial correlation for unordered categorical vars
matched_cov_bal_plot <- all_cov_bal(matched_pop_subset, w_subset, c_unordered_vars = subset(c_subset, select = zip_unordered_cat_var_names),
            ci_appr = "matching", all_cov_names = colnames(c_subset), title = paste("Set of", format(n_rows, scientific = F), "observations"))
ggsave(paste0(dir_proj, "results/covariate_balance/matched_pop_", n_rows, "rows", modifications, ".png"), matched_cov_bal_plot)

# print summary statistics for pseudopopulation counter
summarize_pseudo_counter(matched_pop_subset)

# pseudopopulation, including individual-level covariates (i.e., strata), trimming away unmatched data
matched_obs <- matched_pop_subset$pseudo_pop$row_index[matched_pop_subset$pseudo_pop$counter > 0]
matched_indiv_vars <- indiv_vars_subset[matched_obs] # To Do: consider including ffs_entry_year/ADRD_year in GPS or outcome model
matched_offset <- offset_subset[matched_obs]
matched_data <- cbind(matched_pop_subset$pseudo_pop[matched_pop_subset$pseudo_pop$counter > 0, ], matched_indiv_vars, matched_offset) # to do: check that rows are in same order
matched_data <- as.data.frame(matched_data)

# Examine distribution of exposure and ZIP-level covariates (which were used to match) in pseudopopulation
summary(matched_data$w)
explore_zip_covs(matched_data)



##### IPTW by GPS using CausalGPS package #####

# create log file to see internal processes of CausalGPS
set_logger(logger_file_path = paste0(dir_proj, "code/analysis/CausalGPS_logs/CausalGPS_", Sys.Date(), "_weight_", modifications, n_rows, "rows_", n_cores, "cores_", n_gb, "gb.log"),
           logger_level = "DEBUG")

# GPS weighting on ZIP-level covariates
set.seed(200)
weighted_pop_subset <- generate_pseudo_pop(Y_subset,
                                          w_subset,
                                          c_subset,
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
                                          trim_quantiles = c(0.01,0.99), # c(0.05, 0.95) or c(0.01, 0.99)
                                          max_attempt = 10,
                                          matching_fun = "matching_l1",
                                          delta_n = 0.2, # std dev of pm2.5 is 2.87, so I'll set delta_n = 0.2? parameters may depend on if state-level or national
                                          scale = 1)
saveRDS(weighted_pop_subset, file = paste0(dir_proj, "data/pseudopops/weighted_pop_", n_rows, "rows", modifications, ".rds"))

# check ZIP-level covariate balance
# i.e., absolute correlation for quantitative covariates, polyserial correlation for ordered categorical variables, mean absolute point-biserial correlation for unordered categorical vars
weighted_cov_bal_plot <- all_cov_bal(weighted_pop_subset, w_subset, c_unordered_vars = subset(c_subset, select = zip_unordered_cat_var_names),
            "weighting", all_cov_names = colnames(c_subset), title = paste("Set of", format(n_rows, scientific = F), "observations"))
ggsave(paste0(dir_proj, "results/covariate_balance/weighted_pop_", n_rows, "rows", modifications, ".png"), weighted_cov_bal_plot)

# print summary statistics for pseudopopulation weights
summarize_pseudo_weights(weighted_pop_subset)

# pseudopopulation, including individual-level covariates (i.e., strata), trimming away unmatched data
weighted_obs <- weighted_pop_subset$pseudo_pop$row_index
weighted_indiv_vars <- indiv_vars_subset[weighted_obs] # To Do: consider including ffs_entry_year/ADRD_year in GPS or outcome model
weighted_offset <- offset_subset[weighted_obs]
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
  ipw <- weighted_pop_subset$pseudo_pop$ipw
  capped_weighted_pop_subset <- copy(weighted_pop_subset)
  capped_weighted_pop_subset$pseudo_pop$ipw <- ifelse(ipw > cutoff_weight, cutoff_weight, ipw)
  adjusted_corr_obj <- check_covar_balance(capped_weighted_pop_subset$pseudo_pop,
                                           ci_appr = "weighting",
                                           nthread = n_cores - 1,
                                           covar_bl_method = "absolute",
                                           covar_bl_trs = 0.2, # or 0.1
                                           covar_bl_trs_type = "maximal",
                                           optimized_compile = T)
  capped_weighted_pop_subset$adjusted_corr_results <- adjusted_corr_obj$corr_results
}

cat("ESS of capped weighted pseudopopulation:", ess(capped_weighted_pop_subset$pseudo_pop$ipw))

# check ZIP-level covariate balance
# i.e., absolute correlation for quantitative covariates, polyserial correlation for ordered categorical variables, mean absolute point-biserial correlation for unordered categorical vars
capped_weighted_cov_bal_plot <- all_cov_bal(capped_weighted_pop_subset, w_subset, c_unordered_vars = subset(c_subset, select = zip_unordered_cat_var_names),
            "weighting", all_cov_names = colnames(c_subset), title = paste("Set of", format(n_rows, scientific = F), "observations, weights capped at", cutoff_weight))
ggsave(paste0(dir_proj, "results/covariate_balance/capped_weighted_pop_", n_rows, "rows", modifications, ".png"), capped_weighted_cov_bal_plot)

# pseudopopulation, including individual-level covariates (i.e., strata), trimming away unmatched data
capped_weighted_obs <- capped_weighted_pop_subset$pseudo_pop$row_index
capped_weighted_indiv_vars <- indiv_vars_subset[capped_weighted_obs] # To Do: consider including ffs_entry_year/ADRD_year in GPS or outcome model
capped_weighted_offset <- offset_subset[capped_weighted_obs]
capped_weighted_data <- cbind(capped_weighted_pop_subset$pseudo_pop, capped_weighted_indiv_vars, capped_weighted_offset) # to do: check that rows are in same order
capped_weighted_data <- as.data.frame(capped_weighted_data)

# Examine distribution of exposure and ZIP-level covariates (which were used to match) in pseudopopulation
summary(capped_weighted_data$w)
explore_zip_covs(capped_weighted_data)


##### Adjusting by GPS using CausalGPS package #####

# create log file to see internal processes of CausalGPS
set_logger(logger_file_path = paste0(dir_proj, "code/analysis/CausalGPS_logs/CausalGPS_", Sys.Date(), "_adjust_", modifications, n_rows, "rows_", n_cores, "cores_", n_gb, "gb.log"),
           logger_level = "DEBUG")

# GPS estimation on ZIP-level covariates
set.seed(200)
estimating_gps <- estimate_gps(Y_subset,
                              w_subset,
                              c_subset,
                              pred_model = "sl",
                              gps_model = "parametric",
                              sl_lib = c("m_xgboost"),
                              params = list(xgb_nrounds = seq(10, 50),
                                            xgb_eta = seq(0.1, 0.4, 0.01)),
                              nthread = n_cores - 1,
                              internal_use = FALSE)
estimated_gps <- estimating_gps$dataset
estimated_gps$counter <- NULL # irrelevant to adjusting approach
estimated_gps$row_index <- NULL # irrelevant to adjusting approach
saveRDS(estimated_gps, file = paste0(dir_proj, "data/pseudopops/estimated_gps_", n_rows, "rows", modifications, ".rds"))


# check ZIP-level covariate balance (note: this is the original, unchanged population)
# i.e., absolute correlation for quantitative covariates, polyserial correlation for ordered categorical variables, mean absolute point-biserial correlation for unordered categorical vars
abs_cor_orig <- rep(NA, ncol(c_subset))
names(abs_cor_orig) <- colnames(c_subset)
for (ordered_var in zip_quant_var_names){
  abs_cor_orig[ordered_var] <- abs(cor(w_subset, c_subset[[ordered_var]]))
}
for (unordered_var in zip_unordered_cat_var_names){
  abs_cor_orig[unordered_var] <- abs(cor_unordered_var(w_subset, c_subset[[unordered_var]]))
}
abs_cor_orig <- data.frame(cov = names(abs_cor_orig), abs_cor = abs_cor_orig)
ggplot(abs_cor_orig, aes(x = cov, y = abs_cor)) +
  geom_point() +
  geom_line() +
  labs(title = paste("Set of", format(n_rows, scientific = F), "observations"), x = "Covariate", y = "Absolute Correlation") +
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5))


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