#TODO: Convert it to project's convention.

rm(list = ls())
gc()

## 0. Setup --------------------------------------------------------------------
# devtools::install_github("fasrc/CausalGPS", ref="develop")
library(data.table)
library(fst)
library(CausalGPS)
library(mgcv)
library(ggplot2)
library(tidyr)
library(memoise)
library(parallel)

# TODO: Needs a better way to set the working dir.
setwd("/n/dominici_nsaph_l3/Lab/projects/nkhoshnevis_pm_no2_o3-adrd_hosp-medicare-causalgps/pm_no2_o3-adrd_hospitalization-medicare-causalgps/analysis")
readRenviron(".Renviron")
dir_data <- Sys.getenv("ADRD_CAUSALGPS_DATA")
dir_code <- Sys.getenv("ADRD_CAUSALGPS_CODE")

sub_proj <- "sp_230202_stratified"

interim_sp_results <- file.path(getwd(), sub_proj)
if (!dir.exists(interim_sp_results)) {dir.create(interim_sp_results)}

sp_cache <- file.path(interim_sp_results, "cache_files")
if (!dir.exists(sp_cache)) {dir.create(sp_cache)}

sp_log <- file.path(interim_sp_results, "log_files")
if (!dir.exists(sp_log)) {dir.create(sp_log)}

sp_output <- file.path(interim_sp_results, "output_files")
if (!dir.exists(sp_output)) {dir.create(sp_output)}

source(paste0(dir_code, "analysis/helper_functions.R"))

# setup cache on disk
cdb <- cachem::cache_disk(sp_cache)

## 1. Controlling parameters ---------------------------------------------------

n_cores <- 48 # 48 is max of fasse partition, 64 js max of fasse_bigmem partition
n_gb <- 184 # 184 is max of fasse partition, 499 is max of fasse_bigmem partition
# total_n_rows <- nrow(ADRD_agg_lagged)
n_attempts <- 50
exposure_name <- "pm25"
outcome_name <- "n_hosp"
trim_quantiles <- c(0.01, 0.99)

## 2. Load data ----------------------------------------------------------------

load_data_tmp <- function(data_path){

  # read in full data; 34,763,397 rows
  data <- read_fst(data_path,
                   as.data.table = TRUE)

  setnames(data,
           old = c("pct_blk", "pct_owner_occ"),
           new = c("prop_blk", "prop_owner_occ"))

  data[, `:=`(zip = as.factor(zip),
                         year = as.factor(year),
                         cohort = as.factor(cohort),
                         age_grp = as.factor(age_grp),
                         sex = as.factor(sex),
                         race = as.factor(race),
                         dual = as.factor(dual))]

  return(data)
}

m_load_data_tmp <- memoise(load_data_tmp, cache = cdb)

# # read in full data; 34,763,397 rows
ADRD_agg_lagged <- m_load_data_tmp(file.path(dir_data,
                                             "analysis/ADRD_complete_tv.fst"))

# parameters for this computing job
n_total_attempts <- n_attempts # user can set this to a number larger than n_attempts if some attempts with different seeds have already been tried
modifications <- paste0("gps_by_zip_year_", n_attempts, "attempts") # to be used in names of output files, to record how you're tuning the models


## 3. Preprocess data ----------------------------------------------------------

prep_data <- function(data,
                      fields,
                      zip_expos_names,
                      other_expos_names,
                      exposure_name,
                      outcome_name,
                      zip_unordered_cat_var_names,
                      indiv_unordered_cat_var_names,
                      subset_select,
                      unique_select,
                      trim_quantiles
                      ){

  #other_expos_names <- zip_expos_names[zip_expos_names != exposure_name]
  data_subset <- subset(data,
                        select = fields)

  for (var in c(zip_unordered_cat_var_names, indiv_unordered_cat_var_names)){
    data[[var]] <- as.factor(data[[var]])
  }

  colnames(data_subset)[colnames(data_subset) == exposure_name] <- "w"
  zip_year_data <- subset(data_subset,
                          select = subset_select)

  zip_year_data <- unique(zip_year_data, by = unique_select)
  trim_d_q <- quantile(zip_year_data$w, trim_quantiles)
  zip_year_rows_within_range <- zip_year_data$w >= trim_d_q[1] & zip_year_data$w <= trim_d_q[2]
  #zip_year_trimmed <- zip_year_data[zip_year_rows_within_range, ]
  #n_zip_year_rows <- nrow(zip_year_data) # 486,793

  data_within_range <- data_subset$w >= trim_d_q[1] & data_subset$w <= trim_d_q[2]
  data_trimmed <- data_subset[data_within_range, ]
  #n_rows <- nrow(ADRD_agg_lagged_trimmed_1_99) # 34,141,155 for PM1.5; 34,090,022 for NO2; 33,411,373 for ozone
  return(list(data_trimmed, zip_year_data))
}

m_prep_data <- memoise(prep_data,  cache = cdb)

other_expos_names <- zip_expos_names[zip_expos_names != exposure_name]

fields = c(exposure_name,
           outcome_name,
           other_expos_names,
           zip_var_names,
           "zip",
           indiv_var_names,
           offset_var_names)

subset_select <-  c("zip", "year", "w",
                    other_expos_names, zip_var_names)
unique_select <- c("zip", "year")

processed_data<- m_prep_data(data = ADRD_agg_lagged,
                             fields = fields,
                             zip_expos_names = zip_expos_names,
                             other_expos_names = other_expos_names,
                             exposure_name = exposure_name,
                             outcome_name = outcome_name,
                             zip_unordered_cat_var_names = zip_unordered_cat_var_names,
                             indiv_unordered_cat_var_names = indiv_unordered_cat_var_names,
                             subset_select = subset_select,
                             unique_select = unique_select,
                             trim_quantiles = trim_quantiles)

ADRD_agg_lagged_trimmed_1_99 <- processed_data[[1]]
zip_year_data <- processed_data[[2]]

n_rows <- nrow(ADRD_agg_lagged_trimmed_1_99) # 34,141,155 for PM1.5; 34,090,022 for NO2; 33,411,373 for ozone


## 4. GPS Weighting ------------------------------------------------------------

# set up data.frame to check covariate balance for each GPS modeling attempt
### to do: see if zip can be included or if need more memory or something
vars_for_cov_bal <- c(other_expos_names, zip_var_names, indiv_var_names)
cov_bal_weighting <- expand.grid(1:n_attempts,
                                 vars_for_cov_bal,
                                 c("Weighted", "Unweighted"),
                                 100)

# Absolute_Correlation column will be updated
colnames(cov_bal_weighting) <- c("Attempt",
                                 "Covariate",
                                 "Dataset",
                                 "Absolute_Correlation")

cov_bal_weighting <- as.data.table(cov_bal_weighting)

# set up list to store each GPS modeling attempt
gps_for_weighting_list <- vector("list", n_attempts)

# create log file to see internal processes of CausalGPS
set_logger(logger_file_path = file.path(sp_log, "CausalGPS.log "),
           logger_level = "TRACE")

for (i in 1:n_attempts){

  # estimate GPS
  set.seed(i*100)
  temp_zip_year_with_gps <- estimate_gps(Y = 0, # fake Y variable since our outcomes are not at the zip-year level; not used in estimate_gps
                                         w = zip_year_data$w,
                                         c = subset(zip_year_data, select = c("year", other_expos_names, zip_var_names)),
                                         gps_model = "parametric", # i.e., w=f(x)+epsilon, f(x) estimated by xgboost and epsilon is normal
                                         internal_use = TRUE,
                                         params = list(xgb_nrounds = seq(10, 50),
                                                       xgb_eta = seq(0.1, 0.4, 0.01)),
                                         sl_lib = c("m_xgboost"),
                                         nthread = n_cores)

  temp_zip_year_with_gps <- temp_zip_year_with_gps$dataset
  temp_zip_year_with_gps$zip <- zip_year_data$zip

  # stabilize GPS using marginal probability of exposure (modeled normally)
  # and cap extreme weights at 10
  marginal_expos_prob <- dnorm(zip_year_data$w,
                               mean = mean(zip_year_data$w),
                               sd = sd(zip_year_data$w))
  temp_zip_year_with_gps$stabilized_ipw <- marginal_expos_prob /
                                           temp_zip_year_with_gps$gps ## check estimate_gps
  temp_zip_year_with_gps$capped_stabilized_ipw <- ifelse(temp_zip_year_with_gps$stabilized_ipw > 10,
                                                         10,
                                                         temp_zip_year_with_gps$stabilized_ipw)
  gps_for_weighting_list[[i]] <- temp_zip_year_with_gps

  ### to do: compare results with create_weighting()

  # merge with patient data
  temp_weighted_pseudopop <- merge(ADRD_agg_lagged_trimmed_1_99,
                                   subset(temp_zip_year_with_gps,
                                          select = c("zip",
                                                     "year",
                                                     "capped_stabilized_ipw")),
                                   by = c("zip", "year"))

  # calculate covariate balance: mean absolute point-biserial correlation for unordered categorical vars
  for (unordered_var in c(zip_unordered_cat_var_names)){
    cov_bal_weighting[Attempt == i &
                      Covariate == unordered_var &
                      Dataset == "Weighted",
                      Absolute_Correlation := weighted_cor_unordered_var(temp_weighted_pseudopop$w,
                                                                         temp_weighted_pseudopop[[unordered_var]],
                                                                         temp_weighted_pseudopop$capped_stabilized_ipw)]
    cov_bal_weighting[Attempt == i &
                      Covariate == unordered_var &
                      Dataset == "Unweighted",
                      Absolute_Correlation := cor_unordered_var(temp_weighted_pseudopop$w,
                                                                temp_weighted_pseudopop[[unordered_var]])]
  }

  ## to do: age_grp is ordered

  # calculate covariate balance: absolute correlation for quantitative covariates
  for (quant_var in c(other_expos_names, zip_quant_var_names, indiv_quant_var_names)){
    cov_bal_weighting[Attempt == i &
                      Covariate == quant_var &
                      Dataset == "Weighted",
                      Absolute_Correlation := weightedCorr(temp_weighted_pseudopop$w,
                                                           temp_weighted_pseudopop[[quant_var]],
                                                           method = "Pearson",
                                                           weights = temp_weighted_pseudopop$capped_stabilized_ipw)]
    cov_bal_weighting[Attempt == i &
                      Covariate == quant_var &
                      Dataset == "Unweighted",
                      Absolute_Correlation := cor(temp_weighted_pseudopop$w,
                                                  temp_weighted_pseudopop[[quant_var]],
                                                  method = "pearson")]
    cov_bal_weighting[Attempt == i &
                      Covariate == quant_var,
                      Absolute_Correlation := abs(Absolute_Correlation)]
  }
}


rm(temp_zip_year_with_gps, temp_weighted_pseudopop)

# find GPS model with best covariate balance
# We do not compute covariate balance for sex, race, dual.
# Also there are missing values.
cov_bal_weighting_2 <- cov_bal_weighting[
                        !(cov_bal_weighting$Absolute_Correlation > 1 |
                        is.na(cov_bal_weighting$Absolute_Correlation))]

cov_bal_summary <- cov_bal_weighting_2[Dataset == "Weighted",
                                     .(max_abs_cor = max(Absolute_Correlation,
                                                         na.rm = TRUE)),
                                     by = Attempt]
best_attempt <- cov_bal_summary$Attempt[which.min(cov_bal_summary$max_abs_cor)]
best_cov_bal <- cov_bal_weighting_2[Attempt == best_attempt]

# plot covariate balance
weighted_cov_bal_plot <- ggplot(best_cov_bal,
                                aes(x = Covariate,
                                    y = Absolute_Correlation,
                                    color = Dataset, group = Dataset)) +
                         geom_point() +
                         geom_line() +
  ylab(paste("Absolute Correlation with", exposure_name)) +
  ggtitle(paste0(format(n_rows, scientific = F, big.mark = ','),
                 " units of analysis (Attempt #",
                 best_attempt, " of ",
                 n_total_attempts, ")")) +
  theme(axis.text.x = element_text(angle = 90),
        plot.title = element_text(hjust = 0.5))

ggsave(paste0(sp_output, "/weighted_pop_",
                 n_rows, "_rows_", modifications, ".png"),
       weighted_cov_bal_plot, width = 10, height = 10, dpi = 300)

# print summary statistics for pseudopopulation weights
best_weighted_pseudopop <- merge(ADRD_agg_lagged_trimmed_1_99,
                                 subset(gps_for_weighting_list[[best_attempt]],
                                        select = c("zip",
                                                   "year",
                                                   "capped_stabilized_ipw")),
                                 by = c("zip", "year"))
ess(best_weighted_pseudopop$capped_stabilized_ipw) # 7,889,768

# model outcome from GPS-weighted pseudo-population
formula_expos_only <- as.formula(paste(outcome_name, "~", paste(c("w", strata_vars), collapse = "+", sep = "")))
formula_expos_only_smooth <- as.formula(paste(outcome_name, "~", paste(c("s(w, bs = 'ts')", strata_vars), collapse = "+", sep = "")))

# parametric model (Poisson regression)

cl <- parallel::makeCluster(n_cores, type = "PSOCK")

bam_exposure_only_capped_weighted <- bam(formula_expos_only,
                                         data = best_weighted_pseudopop,
                                         offset = log(n_persons * n_years),
                                         family = poisson(link = "log"),
                                         weights = capped_stabilized_ipw,
                                         samfrac = 0.05,
                                         chunk.size = 5000,
                                         control = gam.control(trace = TRUE),
                                         nthreads = n_cores,
                                         cluster = cl)

parallel::stopCluster(cl)


summary(bam_exposure_only_capped_weighted)
saveRDS(summary(bam_exposure_only_capped_weighted),
        file = paste0(sp_output,
                      "/bam_capped_weighted_exposure_only_",
                      n_rows,
                      "rows_",
                      modifications,
                      ".rds"))

# semi-parametric model (thin-plate spline)
cl <- parallel::makeCluster(n_cores, type = "PSOCK")
bam_smooth_exposure_only_capped_weighted <- bam(formula_expos_only_smooth,
                                                data = best_weighted_pseudopop,
                                                offset = log(n_persons * n_years),
                                                family = poisson(link = "log"),
                                                weights = capped_stabilized_ipw,
                                                samfrac = 0.05,
                                                chunk.size = 5000,
                                                control = gam.control(trace = TRUE),
                                                nthreads = n_cores,
                                                cluster = cl)
parallel::stopCluster(cl)

png(paste0(sp_output, "/semiparametric_results_ERFs_bam_smooth_exposure_only_capped_weighted_",
           n_rows, "rows_", modifications, ".png"))
plot(bam_smooth_exposure_only_capped_weighted,
     main = paste0("GPS-Weighted, Capped at 10, Smoothed Poisson regression,\nexposure only (", exposure_name, ")"))
dev.off()

saveRDS(bam_smooth_exposure_only_capped_weighted,
        file = paste0(sp_output,
                      "/semiparametric_results_spline_objects_bam_smooth_exposure_only_capped_weighted_",
                      n_rows, "rows_", modifications, ".rds"))

run_the_rest_of_code = FALSE

if (run_the_rest_of_code){

## 5. GPS Matching -------------------------------------------------------------

##### THE CODE BELOW IS UNFINISHED #####

# get columns from full data that are useful for matching
ADRD_agg_for_matching <- copy(ADRD_agg_lagged_trimmed_1_99)
ADRD_agg_for_matching[, stratum := .GRP, by = strata_vars]
setorder(ADRD_agg_for_matching, stratum) # to do: consider if this line is necessary
ADRD_agg_for_matching <- ADRD_agg_for_matching[, .(zip, year, stratum)]

# set up dataframe to check covariate balance for each GPS modeling attempt
### to do: see if zip can be included or if need more memory or something
vars_for_cov_bal <- c(other_expos_names, zip_var_names, indiv_var_names)
cov_bal_matching <- expand.grid(1:n_attempts,
                                vars_for_cov_bal,
                                c("Matched", "Unmatched"),
                                100)
colnames(cov_bal_matching) <- c("Attempt", "Covariate", "Dataset", "Absolute_Correlation") # Absolute_Correlation column will be updated
cov_bal_matching <- as.data.table(cov_bal_matching)

# set up list to store each GPS modeling attempt
gps_for_matching_list <- vector("list", n_attempts)

# use same exposure bin sequence for all strata's matching
# bin_seq_by_quantile <- quantile(zip_year_data$w, 0:100/100)
# matching_caliper <- mean(diff(bin_seq_by_quantile))
matching_caliper <- 0.6

for (i in 1:n_attempts){
  # create log file to see internal processes of CausalGPS
  set_logger(logger_file_path = paste0(dir_code, "analysis/CausalGPS_logs/CausalGPS_", Sys.Date(), "_estimate_gps_for_matching_", modifications, "_", n_zip_year_rows, "rows_", n_cores, "cores_", n_gb, "gb.log"),
             logger_level = "TRACE")

  # estimate GPS
  set.seed(i*100)
  temp_zip_year_with_gps_obj <- estimate_gps(Y = 0, # fake Y variable since our outcomes are not at the zip-year level; not used in estimate_gps
                                         w = zip_year_data$w,
                                         c = subset(zip_year_data, select = c("year", other_expos_names, zip_var_names)),
                                         gps_model = "parametric", # i.e., w=f(x)+epsilon, f(x) estimated by xgboost and epsilon is normal
                                         internal_use = T,
                                         params = list(xgb_nrounds = seq(10, 50),
                                                       xgb_eta = seq(0.1, 0.4, 0.01)),
                                         sl_lib = c("m_xgboost"),
                                         nthread = n_cores)

  # create a temporary dataset storing all of CausalGPS's internal parameters, to be expanded from ZIP-years to units of analysis (merging with strata by ZIP, year)
  temp_zip_year_with_gps_dataset_plus_params <- temp_zip_year_with_gps_obj$dataset
  temp_zip_year_with_gps_dataset_plus_params$e_gps_pred <- temp_zip_year_with_gps_obj$e_gps_pred
  temp_zip_year_with_gps_dataset_plus_params$w_resid <- temp_zip_year_with_gps_obj$w_resid

  # apply estimated GPS value to all strata within each ZIP/year
  temp_zip_year_with_gps_dataset_plus_params$zip <- zip_year_data$zip
  temp_zip_year_with_gps_dataset_plus_params <- merge(ADRD_agg_for_matching, temp_zip_year_with_gps_dataset_plus_params,
                                                      by = c("zip", "year"))
  setorder(temp_zip_year_with_gps_dataset_plus_params, stratum)
  temp_zip_year_with_gps_dataset_plus_params[, row_index := 1:.N, by = stratum]
  # temp_zip_year_with_gps_dataset_plus_params$row_index = 1:n_rows # old code
  temp_zip_year_with_gps_dataset_plus_params$zip <- NULL # ZIP is the only free variable for matching, so remove it from the dataset

  # match within strata
  strata_list <- split(temp_zip_year_with_gps_dataset_plus_params, temp_zip_year_with_gps_dataset_plus_params$stratum)
  set_logger(logger_file_path = paste0(dir_code, "analysis/CausalGPS_logs/CausalGPS_", Sys.Date(), "_matching_by_stratum", modifications, "_", n_rows, "rows_", n_cores, "cores_", n_gb, "gb.log"),
             logger_level = "TRACE")

  match_within_stratum <- function(dataset_plus_params){
    # make cgps_gps object from input ("dataset_plus_params")
    dataset_as_cgps_gps <- list()
    class(dataset_as_cgps_gps) <- "cgps_gps"
    dataset_as_cgps_gps$dataset <- subset(dataset_plus_params,
                                          select = c("Y", "w", "year", zip_var_names, # note that Y is a fake variable with value 0; not used in matching
                                                     "gps", "counter_weight", "row_index"))
    # dataset_as_cgps_gps$used_params <- used_params # old code
    dataset_as_cgps_gps$e_gps_pred <- dataset_plus_params$e_gps_pred
    dataset_as_cgps_gps$w_resid <- dataset_plus_params$w_resid

    dataset_as_cgps_gps$e_gps_std_pred <- temp_zip_year_with_gps_obj$e_gps_std_pred
    dataset_as_cgps_gps$gps_mx <- temp_zip_year_with_gps_obj$gps_mx
    dataset_as_cgps_gps$w_mx <- temp_zip_year_with_gps_obj$w_mx

    # the following two approaches to matching throw errors. how to fix?
    matched_pop <- compile_pseudo_pop(data_obj = dataset_as_cgps_gps,
                                      ci_appr = "matching",
                                      gps_model = "parametric",
                                      bin_seq = NULL,
                                      nthread = n_cores,
                                      optimized_compile = T,
                                      matching_fun = "matching_l1",
                                      covar_bl_method = "absolute",
                                      covar_bl_trs = 0.1,
                                      covar_bl_trs_type = "maximal",
                                      delta_n = matching_caliper,
                                      scale = 1) # max_attempt is not a parameter

    matched_pop2 <- CausalGPS:::create_matching(dataset = dataset_as_cgps_gps,
                                                bin_seq = NULL,
                                                gps_model = "parametric",
                                                nthread = n_cores,
                                                optimized_compile = T,
                                                matching_fun = "matching_l1",
                                                delta_n = matching_caliper)

    return()
  }
  # temp_matched_pseudopops <- lapply(strata_list, match_within_stratum)

}

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
saveRDS(matched_pop_subset, file = paste0(dir_data, "pseudopops/matched_pop_", n_rows, "rows", modifications, ".rds"))

# check ZIP-level covariate balance
# i.e., absolute correlation for quantitative covariates, polyserial correlation for ordered categorical variables, mean absolute point-biserial correlation for unordered categorical vars
matched_cov_bal_plot <- all_cov_bal(matched_pop_subset, w, c_unordered_vars = subset(c, select = zip_unordered_cat_var_names),
            ci_appr = "matching", all_cov_names = colnames(c), title = paste("Set of", format(n_rows, scientific = F), "observations"))
ggsave(paste0(dir_results, "covariate_balance/matched_pop_", n_rows, "rows", modifications, ".png"), matched_cov_bal_plot)

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
set_logger(logger_file_path = paste0(dir_code, "analysis/CausalGPS_logs/CausalGPS_", Sys.Date(), "_weight_", modifications, "_", n_rows, "rows_", n_cores, "cores_", n_gb, "gb.log"),
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
saveRDS(weighted_pop_subset, file = paste0(dir_data, "pseudopops/weighted_pop_", n_rows, "rows", modifications, ".rds"))

# check ZIP-level covariate balance
# i.e., absolute correlation for quantitative covariates, polyserial correlation for ordered categorical variables, mean absolute point-biserial correlation for unordered categorical vars
weighted_cov_bal_plot <- all_cov_bal(weighted_pop_subset, w, c_unordered_vars = subset(c, select = zip_unordered_cat_var_names),
            "weighting", all_cov_names = colnames(c), title = paste("Set of", format(n_rows, scientific = F), "observations"))
ggsave(paste0(dir_results, "covariate_balance/weighted_pop_", n_rows, "rows", modifications, ".png"), weighted_cov_bal_plot)

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
ggsave(paste0(dir_results, "covariate_balance/capped_weighted_pop_", n_rows, "rows", modifications, ".png"), capped_weighted_cov_bal_plot)

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
ggsave(paste0(dir_results, "covariate_balance/weighted_pop_", n_rows, "rows", modifications, ".png"), weighted_cov_bal_plot)

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
ggsave(paste0(dir_results, "covariate_balance/capped_weighted_pop_", n_rows, "rows", modifications, ".png"), capped_weighted_cov_bal_plot)


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


}
