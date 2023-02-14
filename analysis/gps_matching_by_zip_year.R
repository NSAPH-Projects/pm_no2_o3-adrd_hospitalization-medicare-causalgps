dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/code/"
source(paste0(dir_code, "analysis/setup_trimmed_and_zip_year_data.R"))


##### GPS Matching #####

# get columns from full data that are useful for matching
ADRD_agg_for_matching <- copy(ADRD_agg_lagged_trimmed_1_99)
ADRD_agg_for_matching[, stratum := .GRP, by = strata_vars]
setorder(ADRD_agg_for_matching, stratum) # to do: consider if this line is necessary
ADRD_agg_for_matching <- ADRD_agg_for_matching[, .(zip, year, stratum)]

# set up dataframe to check covariate balance for each GPS modeling attempt
### to do: see if zip can be included or if need more memory or something
vars_for_cov_bal <- c(other_expos_names, zip_quant_var_names, levels(zip_year_data[["region"]]))
cov_bal_matching <- expand.grid(1:n_attempts,
                                vars_for_cov_bal,
                                c("Matched", "Unmatched"),
                                100,
                                100)
colnames(cov_bal_matching) <- c("Attempt", "Covariate", "Dataset", "Correlation", "Absolute_Correlation") # Correlation and Absolute_Correlation columns will be updated
cov_bal_matching <- as.data.table(cov_bal_matching)

# use same exposure bin sequence for all strata's matching
matching_caliper <- 0.6 ## to do: play with this
# bin_seq_by_quantile <- quantile(zip_year_data$w, 0:100/100)
# matching_caliper <- mean(diff(bin_seq_by_quantile))

# function to return a pseudopopulation using data from 1 stratum
# as a data.table with variable "counter_weight" denoting number of times each observation is matched
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

# function to match within every stratum and append all those resulting pseudopopulations into a single "full" pseudopopulation
# for one attempt at GPS estimation (a unique attempt_number corresponds to a unique seed)
# can return the covariate balance (as a data.table) OR the "full" pseudopopulation (as a list of pseudopopulations)
get_matched_pseudopop <- function(attempt_number,
                                  cov_bal_matrix,
                                  return_cov_bal = T,
                                  return_pseudopop_list = F){
  
  # create log file to see internal processes of CausalGPS
  set_logger(logger_file_path = paste0(dir_code, "analysis/CausalGPS_logs/CausalGPS_", Sys.Date(), "_estimateGpsForMatching_AttemptNumber", attempt_number, "_", modifications, "_", n_zip_year_rows, "rows_", n_cores, "cores_", n_gb, "gb.log"),
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
  temp_zip_year_with_gps_dataset_plus_params <- merge(ADRD_agg_for_matching, temp_zip_year_with_gps_dataset_plus_params,
                                                      by = c("zip", "year"))
  setorder(temp_zip_year_with_gps_dataset_plus_params, stratum)
  temp_zip_year_with_gps_dataset_plus_params[, row_index := 1:.N, by = stratum]
  # temp_zip_year_with_gps_dataset_plus_params$row_index = 1:n_rows # old code
  temp_zip_year_with_gps_dataset_plus_params$zip <- NULL # ZIP is the only free variable for matching, so remove it from the dataset
  
  # match within strata
  strata_list <- split(temp_zip_year_with_gps_dataset_plus_params,
                       temp_zip_year_with_gps_dataset_plus_params$stratum)
  set_logger(logger_file_path = paste0(dir_code, "analysis/CausalGPS_logs/CausalGPS_", Sys.Date(), "_matching_by_stratum", modifications, "_", n_rows, "rows_", n_cores, "cores_", n_gb, "gb.log"),
             logger_level = "TRACE")
  temp_matched_pseudopop_list <- lapply(strata_list,
                                        match_within_stratum,
                                        e_gps_std_pred = temp_zip_year_with_gps_obj$e_gps_std_pred,
                                        gps_mx = temp_zip_year_with_gps_obj$gps_mx,
                                        w_mx = temp_zip_year_with_gps_obj$w_mx)
  
  if (return_cov_bal){
    temp_matched_pseudopop <- rbindlist(temp_matched_pseudopop_list)
    
    for (unordered_var in c(zip_unordered_cat_var_names)){
      for (level in levels(zip_year_data[[unordered_var]])){
        cov_bal_matrix[Attempt == i & Covariate == level & Dataset == "Matched", Correlation := weightedCorr(temp_matched_pseudopop$w, temp_matched_pseudopop[[unordered_var]] == level, method = "Pearson", weights = temp_matched_pseudopop$counter_weight)]
        cov_bal_matrix[Attempt == i & Covariate == level & Dataset == "Unmatched", Correlation := cor(temp_matched_pseudopop$w, temp_matched_pseudopop[[unordered_var]] == level, method = "pearson")]
      }
    }
    
    for (quant_var in c(other_expos_names, zip_quant_var_names)){
      cov_bal_matrix[Attempt == i & Covariate == quant_var & Dataset == "Matched", Correlation := weightedCorr(temp_matched_pseudopop$w, temp_matched_pseudopop[[quant_var]], method = "Pearson", weights = temp_matched_pseudopop$counter_weight)]
      cov_bal_matrix[Attempt == i & Covariate == quant_var & Dataset == "Unmatched", Correlation := cor(temp_matched_pseudopop$w, temp_matched_pseudopop[[quant_var]], method = "pearson")]
    }
    
    return(cov_bal_matrix)
  }
  
  else if (return_pseudopop_list) return(temp_matched_pseudopop_list)
  else stop("User must specify whether to return covariate balance or pseudopopulations")
}

for (i in 1:n_attempts){
  cov_bal_matching <- get_matched_pseudopop(attempt_number = i,
                                            cov_bal_matrix = cov_bal_matching,
                                            return_cov_bal = T,
                                            return_pseudopop_list = F)
}

# identify GPS model with best covariate balance
cov_bal_matching[, Absolute_Correlation := abs(Correlation)]
cov_bal_matching_summary <- cov_bal_matching[Dataset == "Matched", .(max_abs_cor = max(Absolute_Correlation)), by = Attempt]
best_matching_attempt <- cov_bal_matching_summary$Attempt[which.min(cov_bal_matching_summary$max_abs_cor)]
best_matching_cov_bal <- cov_bal_matching[Attempt == best_matching_attempt]

# plot covariate balance
matched_cov_bal_plot <- ggplot(best_matching_cov_bal, aes(x = Covariate, y = Absolute_Correlation, color = Dataset, group = Dataset)) +
  geom_point() +
  geom_line() +
  ylab(paste("Absolute Correlation with", exposure_name)) +
  ggtitle(paste0(format(n_rows, scientific = F, big.mark = ','), " units of analysis (Attempt #", best_matching_attempt, " of ", n_total_attempts, ")")) +
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5))

ggsave(paste0(dir_results, "covariate_balance/matched_pop_", n_rows, "rows_", modifications, ".png"), matched_cov_bal_plot)

# regenerate GPS model and matched pseudopopulation with best covariate balance
# add "stratum" variable back to pseudopop, so that individual variables can be merged back in, for outcome modeling
best_matched_pseudopop_list <- get_matched_pseudopop(attempt_number = best_matching_attempt,
                                                     cov_bal_matrix = cov_bal_matching,
                                                     return_cov_bal = F,
                                                     return_pseudopop_list = T)
temp_matched_pseudopop_list <- lapply(1:length(temp_matched_pseudopop_list), function(i) temp_matched_pseudopop_list[[i]][, stratum := i])
best_matched_pseudopop <- rbindlist(temp_matched_pseudopop_list)

# print summary statistics for pseudopopulation weights
cat("ESS:", ess(best_matched_pseudopop$counter_weight)) # to do: if ESS is small, investigate which observation(s) are being matched so many times and if increasing? or changing caliper helps
cat("Number of observations matched:", sum(best_matched_pseudopop$counter_weight > 0))
cat("Distribution of number of matches per observations:")
summary(best_matched_pseudopop$counter_weight)

# model outcome from GPS-matched pseudo-population
formula_expos_only <- as.formula(paste(outcome_name, "~", paste(c("w", strata_vars), collapse = "+", sep = "")))
formula_expos_only_smooth <- as.formula(paste(outcome_name, "~", paste(c("s(w, bs = 'ts')", strata_vars), collapse = "+", sep = "")))

# parametric model (Poisson regression)
bam_exposure_only_matched <- bam(formula_expos_only,
                                         data = best_matched_pseudopop,
                                         offset = log(n_persons * n_years),
                                         family = poisson(link = "log"),
                                         weights = counter_weight,
                                         samfrac = 0.05,
                                         chunk.size = 5000,
                                         control = gam.control(trace = TRUE))
summary(bam_exposure_only_matched)
saveRDS(summary(bam_exposure_only_matched), file = paste0(dir_results, "parametric_results/bam_matched_exposure_only_", n_rows, "rows_", modifications, ".rds"))

# semi-parametric model (thin-plate spline)
bam_smooth_exposure_only_matched <- bam(formula_expos_only_smooth,
                                                data = best_matched_pseudopop,
                                                offset = log(n_persons * n_years),
                                                family = poisson(link = "log"),
                                                weights = counter_weight,
                                                samfrac = 0.05,
                                                chunk.size = 5000,
                                                control = gam.control(trace = TRUE),
                                                nthreads = n_cores - 1)
png(paste0(dir_results, "semiparametric_results/ERFs/bam_smooth_exposure_only_matched_", n_rows, "rows_", modifications, ".png"))
plot(bam_smooth_exposure_only_matched, main = paste0("GPS-Matched, Smoothed Poisson regression,\nexposure only (", exposure_name, ")"))
dev.off()
saveRDS(bam_smooth_exposure_only_matched, file = paste0(dir_results, "semiparametric_results/spline_objects/bam_smooth_exposure_only_matched_", n_rows, "rows_", modifications, ".rds"))

