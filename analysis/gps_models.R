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

dir_proj <- "~/nsaph_projects/pm_no2_o3-adrd_hosp-medicare-causalgps/"
dir_data <- paste0(dir_proj, "data/")
dir_git <- paste0(dir_proj, "pm_no2_o3-adrd_hospitalization-medicare-causalgps/")

# read in full data
ADRD_agg <- read_fst(paste0(dir_proj, "data/analysis/ADRD_complete_corrected.fst"), as.data.table = TRUE)
ADRD_agg_lagged <- ADRD_agg[ADRD_year - ffs_entry_year >= 2, ] # Approximate first ADRD hospitalization by requiring no ADRD hosps for 2 years

# to do: change name of repo from "pm_no2_o3-adrd_hospitalization-medicare-causalgps" (since name is similar to project folder) to "git"
source(paste0(dir_proj, "pm_no2_o3-adrd_hospitalization-medicare-causalgps/analysis/helper_functions.R"))


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


##### Use (random) subset of data to make code run faster than full data #####

set.seed(100)
n_random_rows <- 100000 # To Do: use a bigger sample
random_rows <- sample(1:nrow(ADRD_agg_lagged_subset), n_random_rows)
Y_subset <- ADRD_agg_lagged_subset[random_rows, n_ADRDhosp]
w_subset <- ADRD_agg_lagged_subset[random_rows, pm25]
c_subset <- as.data.frame(subset(ADRD_agg_lagged_subset[random_rows,], select = zip_var_names))

# Not used in GPS matching, but used in outcome model
indiv_vars_subset <- ADRD_agg_lagged_subset[random_rows, .(sexM, race_cat, any_dual, ADRD_age)] # To Do: consider including ffs_entry_year/ADRD_year in GPS or outcome model
offset_subset <- ADRD_agg_lagged_subset[random_rows, .(person_years = n_persons * n_years)]
data_subset <- cbind(Y_subset, w_subset, indiv_vars_subset, c_subset, offset_subset)

# To Do at later stage: Full data
# Y <- ADRD_agg_lagged$n_ADRDhosp
# w <- ADRD_agg_lagged$pm25
# c <- as.data.frame(subset(ADRD_agg_lagged, select = zip_var_names))


##### Match on GPS using CausalGPS package #####

# create log file to see internal processes of CausalGPS
set_logger(logger_file_path = paste0(dir_git, "analysis/CausalGPS_logs/CausalGPS_09Aug_xgb_1e5rows_32cores_64gb.log"),
           logger_level = "DEBUG")

# if using SL.gam, remove mgcv library and allow custom parameters
# library(SuperLearner)
# m_gam <- function(cts.num = 16, deg.gam = 1, ...) SL.gam(cts.num = cts.num, deg.gam = deg.gam, ...) # to do: explain these params
# detach("package:mgcv", unload=TRUE)

# To Do: consider using larger delta_n
# GPS matching by ZIP-level covariates
set.seed(200)
matched_pop_subset <- generate_pseudo_pop(Y_subset,
                                 w_subset,
                                 c_subset,
                                  ci_appr = "matching",
                                  pred_model = "sl",
                                  gps_model = "parametric",
                                  use_cov_transform = TRUE,
                                  transformers = list("pow2", "pow3", "sqrt", "log_nonneg", "logit_nonneg"),
                                  sl_lib = c("m_xgboost"), # or SL.glm
                                  params = list(xgb_nrounds = c(10, 20, 30, 50)), # comment out if using sl_lib = "SL.glm"
                                  nthread = 31, # 47
                                  covar_bl_method = "absolute",
                                  covar_bl_trs = 0.1,
                                  covar_bl_trs_type = "maximal",
                                 optimized_compile = TRUE,
                                  trim_quantiles = c(0.05,0.95),
                                  max_attempt = 10, # to do: try up to 20
                                  matching_fun = "matching_l1",
                                  delta_n = 0.1, # to do: check if this is good or if should use bigger
                                  scale = 1)

# check ZIP-level covariate balance
# i.e., absolute correlation for quantitative covariates, polyserial correlation for ordered categorical variables, mean absolute point-biserial correlation for unordered categorical vars
all_cov_bal(matched_pop_subset, w_subset, subset(c_subset, select = zip_unordered_cat_var_names),
            "matching", colnames(c_subset), title = paste("Random set of", n_random_rows, "observations"))

# print summary statistics for pseudopopulation counter
summarize_pseudo_weights(matched_pop_subset, "matching")

### old code ###

# check number of matches
counter <- matched_pop_subset$pseudo_pop$counter
# cat("Number of observations included in pseudo-population:", length(unique(matched_pop_subset$pseudo_pop$row_index)))
# cat("Total number of matches:", sum(counter))
# print("Distribution of number of matches per included observation:")
summary(counter)

# check number of observations used in pseudopopulation
cat("Number of observations used in pseudopopulation:", sum(counter > 0))

# check Kish's effective sample size
cat("Kish ESS:", ess(counter))


##### If desired, cap matches ("counter" variable) at 95th percentile #####
# July 27 edit: don't do this for matching

cap_counts = F
if (cap_counts){
  cap_weights(matched_pop_subset, "matching", 15, zip_quant_var_names, zip_unordered_cat_var_names, paste("Random set of", n_random_rows, "observations"))
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
                                          transformers = list("sqrt", "log_nonneg", "logit_nonneg", "pow2", "pow3"), # list("pow2", "pow3")
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

# check ZIP-level covariate balance
# i.e., absolute correlation for quantitative covariates, polyserial correlation for ordered categorical variables, mean absolute point-biserial correlation for unordered categorical vars
all_cov_bal(matched_pop_subset, w_subset, subset(c_subset, select = zip_unordered_cat_var_names),
            "weighting", colnames(c_subset), title = paste("Random set of", n_random_rows, "observations"))

# check number of matches
ipw <- weighted_pop_subset$pseudo_pop$ipw
cat("Number of observations included in pseudo-population:", length(unique(matched_pop_subset$pseudo_pop$row_index)))
cat("Sum of weights:", sum(ipw))
print("Distribution of weights:")
summary(ipw)
quantile(ipw, c(0, 0.25, 0.5, 0.75, 0.95, 0.99, 0.999))

# check Kish's effective sample size
ipw <- weighted_pop_subset$pseudo_pop$ipw
ess(ipw)


##### If desired, cap weights at 95th percentile #####

cap_weights = T

if (cap_weights){
  cutoff_weight <- quantile(weighted_pop_subset$pseudo_pop$ipw, 0.95)
  weighted_pop_subset$pseudo_pop$ipw <- ifelse(weighted_pop_subset$pseudo_pop$ipw > cutoff_weight, cutoff_weight, weighted_pop_subset$pseudo_pop$ipw)
  adjusted_corr_obj <- check_covar_balance(weighted_pop_subset$pseudo_pop,
                                           ci_appr="weighting",
                                           nthread=15,
                                           covar_bl_method = "absolute",
                                           covar_bl_trs = 0.1,
                                           covar_bl_trs_type = "maximal",
                                           optimized_compile=T)
  weighted_pop_subset$adjusted_corr_results <-  adjusted_corr_obj$corr_results
  
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
}


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