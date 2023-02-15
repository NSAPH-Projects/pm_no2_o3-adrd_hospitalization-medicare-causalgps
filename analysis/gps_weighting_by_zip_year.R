# directories for data, code, and results
dir_data <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/data/"
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/code/"
dir_results <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/results/"

source(paste0(dir_code, "analysis/setup_trimmed_and_zip_year_data.R"))

# parameters for this computing job
n_cores <- 8 # 48 is max of fasse partition, 64 js max of fasse_bigmem partition
n_gb <- 64 # 184 is max of fasse partition, 499 is max of fasse_bigmem partition
# total_n_rows <- nrow(ADRD_agg_lagged)
n_attempts <- 10
n_total_attempts <- n_attempts # user can set this to a number larger than n_attempts if some attempts with different seeds have already been tried
modifications <- paste0("gps_by_zip_year_", n_attempts, "attempts") # to be used in names of output files, to record how you're tuning the models


##### GPS Weighting #####

# set up data.table to check covariate balance for each GPS modeling attempt
### to do: see if zip can be included or if need more memory or something
cov_bal_weighting <- create_cov_bal_data.table("weighting", n_attempts)

# set up list to store each GPS modeling attempt
gps_for_weighting_list <- vector("list", n_attempts)

# create log file to see internal processes of CausalGPS
set_logger(logger_file_path = paste0(dir_code, "analysis/CausalGPS_logs/CausalGPS_", Sys.Date(), "_estimate_gps_for_weighting_", modifications, "_", n_zip_year_rows, "rows_", n_cores, "cores_", n_gb, "gb.log"),
           logger_level = "TRACE")

for (i in 1:n_attempts){

  # estimate GPS
  set.seed(i*100)
  temp_zip_year_with_gps <- estimate_gps(Y = 0, # fake Y variable since our outcomes are not at the zip-year level; not used in estimate_gps
                                        w = zip_year_data$w,
                                        c = subset(zip_year_data, select = c("year", other_expos_names, zip_var_names)),
                                        gps_model = "parametric", # i.e., w=f(x)+epsilon, f(x) estimated by xgboost and epsilon is normal
                                        internal_use = T,
                                        params = list(xgb_nrounds = seq(10, 50),
                                                      xgb_eta = seq(0.1, 0.4, 0.01)),
                                        sl_lib = c("m_xgboost"),
                                        nthread = n_cores)
  temp_zip_year_with_gps <- temp_zip_year_with_gps$dataset
  temp_zip_year_with_gps$zip <- zip_year_data$zip
  
  # stabilize GPS using marginal probability of exposure (modeled normally) and cap extreme weights at 10
  marginal_expos_prob <- dnorm(zip_year_data$w,
                               mean = mean(zip_year_data$w),
                               sd = sd(zip_year_data$w))
  temp_zip_year_with_gps$stabilized_ipw <- marginal_expos_prob / temp_zip_year_with_gps$gps ## check estimate_gps
  temp_zip_year_with_gps$capped_stabilized_ipw <- ifelse(temp_zip_year_with_gps$stabilized_ipw > 10, 10, temp_zip_year_with_gps$stabilized_ipw)
  gps_for_weighting_list[[i]] <- temp_zip_year_with_gps
  
  # merge with patient data
  temp_weighted_pseudopop <- merge(ADRD_agg_lagged_trimmed_1_99, subset(temp_zip_year_with_gps,
                                                                   select = c("zip", "year", "capped_stabilized_ipw")),
                              by = c("zip", "year"))

  # calculate correlation between exposure and covariates
  cov_bal_weighting <- calculate_correlations(cov_bal_data.table = cov_bal_weighting,
                                              method = "weighting",
                                              attempt = i,
                                              pseudopop = temp_weighted_pseudopop)
}
rm(temp_zip_year_with_gps, temp_weighted_pseudopop)

# find GPS model(s) with best covariate balance
cov_bal_summary <- summarize_cov_bal(cov_bal_data.table = cov_bal_weighting,
                                        method = "weighting",
                                        save_csv = T)
best_maxAC_attempt <- cov_bal_summary$Attempt[which.min(cov_bal_summary$maxAC)]
best_maxAC_cov_bal <- cov_bal_weighting[Attempt == best_maxAC_attempt]

# plot covariate balance
weighted_cov_bal_plot <- ggplot(best_maxAC_cov_bal, aes(x = Covariate, y = Absolute_Correlation, color = Dataset, group = Dataset)) +
  geom_point() +
  geom_line() +
  ylab(paste("Absolute Correlation with", exposure_name)) +
  ggtitle(paste0(format(n_rows, scientific = F, big.mark = ','), " units of analysis (Attempt #", best_attempt, " of ", n_total_attempts, ")")) +
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5))

ggsave(paste0(dir_results, "covariate_balance/weighted_pop_", n_rows, "rows_", modifications, ".png"), weighted_cov_bal_plot)

# print summary statistics for pseudopopulation weights
best_maxAC_weighted_pseudopop <- merge(ADRD_agg_lagged_trimmed_1_99, subset(gps_for_weighting_list[[best_maxAC_attempt]],
                                                                      select = c("zip", "year", "capped_stabilized_ipw")),
                                 by = c("zip", "year"))
ess(best_maxAC_weighted_pseudopop$capped_stabilized_ipw) # for attempt #121, which is best out of 200, ESS is 9,427,355

# run parametric and semiparametric (thin-plate spline) outcome models
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