##### Setup #####

# devtools::install_github("fasrc/CausalGPS", ref="develop")
library(data.table)
library(fst)
library(CausalGPS)
library(wCorr)
library(mgcv)
library(ggplot2)

# directories for data, code, and results
dir_data <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/data/"
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/code/"
dir_results <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/results/"

# get data and helpful functions
source(paste0(dir_code, "analysis/helper_functions.R"))
zip_year_data <- read.fst(paste0(dir_data, "analysis/pm_zip_year_data_trimmed_1_99.fst"), as.data.table = T)
zip_year_data_with_strata <- read.fst(paste0(dir_data, "analysis/ADRD_agg_lagged_trimmed_1_99.fst"), as.data.table = T)

# parameters for this computing job
n_cores <- 8 # 48 is max of fasse partition, 64 js max of fasse_bigmem partition
n_gb <- 64 # 184 is max of fasse partition, 499 is max of fasse_bigmem partition
n_attempts <- 10
n_total_attempts <- n_attempts # user can set this to a number larger than n_attempts if some attempts with different seeds have already been tried
modifications <- paste0("gps_by_zip_year_", n_attempts, "attempts") # to be used in names of output files, to record how you're tuning the models


##### GPS Matching #####

# get columns from full data that are useful for matching
data_for_matching <- copy(zip_year_data_with_strata)
data_for_matching[, stratum := .GRP, by = strata_vars]
setorder(data_for_matching, stratum) # to do: consider if this line is necessary
data_for_matching <- data_for_matching[, .(zip, year, stratum)]

# set up data.table to check covariate balance for each GPS modeling attempt
### to do: see if zip can be included or if need more memory or something
cov_bal_matching <- create_cov_bal_data.table("matching", n_attempts)

# use same exposure bin sequence for all strata's matching
matching_caliper <- 0.6 ## to do: play with this. goal is to get large enough ESS
# bin_seq_by_quantile <- quantile(zip_year_data$w, 0:100/100)
# matching_caliper <- mean(diff(bin_seq_by_quantile))

# returns a data.table with variable "counter_weight" denoting number of times each observation is matched
match_within_stratum <- function(dataset_plus_params,
                                 e_gps_std_pred,
                                 gps_mx,
                                 w_mx){
  # make cgps_gps object from input ("dataset_plus_params")
  dataset_as_cgps_gps <- list()
  class(dataset_as_cgps_gps) <- "cgps_gps"
  dataset_as_cgps_gps$dataset <- subset(as.data.frame(dataset_plus_params),
                                        select = c("Y", "w", "year", other_expos_names, zip_var_names, # note that Y is a fake variable with value 0; not used in matching
                                                   "gps", "counter_weight", "row_index"))
  # dataset_as_cgps_gps$used_params <- used_params # old code
  dataset_as_cgps_gps$e_gps_pred <- dataset_plus_params$e_gps_pred
  dataset_as_cgps_gps$w_resid <- dataset_plus_params$w_resid
  
  dataset_as_cgps_gps$e_gps_std_pred <- e_gps_std_pred
  dataset_as_cgps_gps$gps_mx <- gps_mx
  dataset_as_cgps_gps$w_mx <- w_mx
  
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
                                    scale = 1) # notice: max_attempt and transformers are not parameters
  
  return(matched_pop)
}

get_matched_pseudopop <- function(attempt_number,
                                  cov_bal_data.table,
                                  return_cov_bal = T,
                                  return_pseudopop_list = F){
  
  # create log file to see internal processes of CausalGPS
  set_logger(logger_file_path = paste0(dir_code, "analysis/CausalGPS_logs/CausalGPS_", Sys.Date(), "_estimateGpsForMatching_AttemptNumber", attempt_number, "_", modifications, "_", nrow(zip_year_data), "rows_", n_cores, "cores_", n_gb, "gb.log"),
             logger_level = "TRACE")
  
  # estimate GPS
  set.seed(attempt_number*100)
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
  temp_zip_year_with_gps_dataset_plus_params <- merge(data_for_matching, temp_zip_year_with_gps_dataset_plus_params,
                                                      by = c("zip", "year"))
  setorder(temp_zip_year_with_gps_dataset_plus_params, stratum)
  temp_zip_year_with_gps_dataset_plus_params[, row_index := 1:.N, by = stratum]
  # temp_zip_year_with_gps_dataset_plus_params$row_index = 1:nrow(zip_year_data_with_strata) # old code
  temp_zip_year_with_gps_dataset_plus_params$zip <- NULL # ZIP is the only free variable for matching, so remove it from the dataset
  
  # match within strata
  strata_list <- split(temp_zip_year_with_gps_dataset_plus_params,
                       temp_zip_year_with_gps_dataset_plus_params$stratum)
  set_logger(logger_file_path = paste0(dir_code, "analysis/CausalGPS_logs/CausalGPS_", Sys.Date(), "_matching_by_stratum", modifications, "_", nrow(zip_year_data_with_strata), "rows_", n_cores, "cores_", n_gb, "gb.log"),
             logger_level = "TRACE")
  temp_matched_pseudopop_list <- lapply(strata_list,
                                        match_within_stratum,
                                        e_gps_std_pred = temp_zip_year_with_gps_obj$e_gps_std_pred,
                                        gps_mx = temp_zip_year_with_gps_obj$gps_mx,
                                        w_mx = temp_zip_year_with_gps_obj$w_mx)
  
  if (return_cov_bal){
    temp_matched_pseudopop <- rbindlist(temp_matched_pseudopop_list)
    cov_bal_data.table <- calculate_correlations(cov_bal_data.table = cov_bal_data.table,
                                                 method = "matching",
                                                 attempt = attempt_number,
                                                 pseudopop = temp_matched_pseudopop)
    return(cov_bal_data.table)
  }
  
  else if (return_pseudopop_list) return(temp_matched_pseudopop_list)
  else stop("User must specify whether to return covariate balance or (list of) pseudopopulations")
}

for (i in 1:n_attempts){
  cov_bal_matching <- get_matched_pseudopop(attempt_number = i,
                                            cov_bal_data.table = cov_bal_matching,
                                            return_cov_bal = T,
                                            return_pseudopop_list = F)
}

# identify GPS model(s) with best covariate balance
cov_bal_summary <- summarize_cov_bal(cov_bal_data.table = cov_bal_matching,
                                     method = "matching",
                                     save_csv = T)
best_maxAC_attempt <- cov_bal_summary$Attempt[which.min(cov_bal_summary$maxAC)]
best_maxAC_cov_bal <- cov_bal_matching[Attempt == best_maxAC_attempt]

# plot covariate balance
matched_cov_bal_plot <- ggplot(best_maxAC_cov_bal, aes(x = Covariate, y = Absolute_Correlation, color = Dataset, group = Dataset)) +
  geom_point() +
  geom_line() +
  ylab(paste("Absolute Correlation with", exposure_name)) +
  ggtitle(paste0(format(nrow(zip_year_data_with_strata), scientific = F, big.mark = ','), " units of analysis (Attempt #", best_maxAC_attempt, " of ", n_total_attempts, ")")) +
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5))

ggsave(paste0(dir_results, "covariate_balance/matched_pop_", nrow(zip_year_data_with_strata), "rows_", modifications, ".png"), matched_cov_bal_plot)

# regenerate GPS model and matched pseudopopulation with best covariate balance
# add "stratum" variable back to pseudopop, so that individual variables can be merged back in, for outcome modeling
best_matched_pseudopop_list <- get_matched_pseudopop(attempt_number = best_maxAC_attempt,
                                                     cov_bal_data.table = cov_bal_matching,
                                                     return_cov_bal = F,
                                                     return_pseudopop_list = T)
best_matched_pseudopop_list <- lapply(1:length(best_matched_pseudopop_list), function(i) best_matched_pseudopop_list[[i]][, stratum := i])
best_matched_pseudopop <- rbindlist(best_matched_pseudopop_list)

# print summary statistics for pseudopopulation weights
cat("ESS:", ess(best_matched_pseudopop$counter_weight)) # to do: if ESS is small, investigate which observation(s) are being matched so many times and if increasing? or changing caliper helps
cat("Number of observations matched:", sum(best_matched_pseudopop$counter_weight > 0))
cat("Proportion of observations matched:", sum(best_matched_pseudopop$counter_weight > 0) / nrow(best_matched_pseudopop))
cat("Distribution of number of matches per observations:")
summary(best_matched_pseudopop$counter_weight)
quantile(best_matched_pseudopop$counter_weight, c(0, 0.25, 0.5, 0.75, 0.95, 0.99, 0.999, 0.9999, 1))
plot(density(best_matched_pseudopop$w,
             weights = best_matched_pseudopop$counter_weight / sum(best_matched_pseudopop$counter_weight)),
     main = "Density of exposure in matched pseudopopulation",
     xlab = exposure_name)
plot(density(best_matched_pseudopop$w),
     main = "Density of exposure in original pseudopopulation",
     xlab = exposure_name)

# run parametric and semiparametric (thin-plate spline) outcome models
weights <- best_matched_pseudopop$counter_weight # note: to use the following functions, need to have "weights" in global environment; to do: improve this
parametric_model_summary <- get_outcome_model_summary(best_matched_pseudopop,
                                                      "matching",
                                                      n_cores,
                                                      "parametric",
                                                      save_results = T)

semiparametric_model_summary <- get_outcome_model_summary(best_matched_pseudopop,
                                                          "matching",
                                                          n_cores,
                                                          "semiparametric",
                                                          save_results = T)