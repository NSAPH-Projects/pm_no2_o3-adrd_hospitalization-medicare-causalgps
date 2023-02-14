dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/code/"
source(paste0(dir_code, "analysis/setup_trimmed_and_zip_year_data.R"))


##### GPS Weighting #####

# set up dataframe to check covariate balance for each GPS modeling attempt
### to do: see if zip can be included or if need more memory or something
vars_for_cov_bal <- c(other_expos_names, zip_quant_var_names, levels(zip_year_data[["region"]]))
cov_bal_weighting <- expand.grid(1:n_attempts,
                                 vars_for_cov_bal,
                                 c("Weighted", "Unweighted"),
                                 100,
                                 100)
colnames(cov_bal_weighting) <- c("Attempt", "Covariate", "Dataset", "Correlation", "Absolute_Correlation") # Correlation and Absolute_Correlation columns will be updated
cov_bal_weighting <- as.data.table(cov_bal_weighting)

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

  for (unordered_var in c(zip_unordered_cat_var_names)){
    for (level in levels(zip_year_data[[unordered_var]])){
      cov_bal_weighting[Attempt == i & Covariate == level & Dataset == "Weighted", Correlation := weightedCorr(temp_weighted_pseudopop$w, temp_weighted_pseudopop[[unordered_var]] == level, method = "pearson", weights = temp_weighted_pseudopop$capped_stabilized_ipw)]
      cov_bal_weighting[Attempt == i & Covariate == level & Dataset == "Unweighted", Correlation := cor(temp_weighted_pseudopop$w, temp_weighted_pseudopop[[unordered_var]] == level, method = "pearson")]
    }
  }
  
  for (quant_var in c(other_expos_names, zip_quant_var_names)){
    cov_bal_weighting[Attempt == i & Covariate == quant_var & Dataset == "Weighted", Correlation := weightedCorr(temp_weighted_pseudopop$w, temp_weighted_pseudopop[[quant_var]], method = "Pearson", weights = temp_weighted_pseudopop$capped_stabilized_ipw)]
    cov_bal_weighting[Attempt == i & Covariate == quant_var & Dataset == "Unweighted", Correlation := cor(temp_weighted_pseudopop$w, temp_weighted_pseudopop[[quant_var]], method = "pearson")]
  }
}
rm(temp_zip_year_with_gps, temp_weighted_pseudopop)

# find GPS model with best covariate balance
cov_bal_weighting[, Absolute_Correlation := abs(Correlation)]
cov_bal_summary <- cov_bal_weighting[Dataset == "Weighted", .(max_abs_cor = max(Absolute_Correlation)), by = Attempt]
best_attempt <- cov_bal_summary$Attempt[which.min(cov_bal_summary$max_abs_cor)]
best_cov_bal <- cov_bal_weighting[Attempt == best_attempt]

# plot covariate balance
weighted_cov_bal_plot <- ggplot(best_cov_bal, aes(x = Covariate, y = Absolute_Correlation, color = Dataset, group = Dataset)) +
  geom_point() +
  geom_line() +
  ylab(paste("Absolute Correlation with", exposure_name)) +
  ggtitle(paste0(format(n_rows, scientific = F, big.mark = ','), " units of analysis (Attempt #", best_attempt, " of ", n_total_attempts, ")")) +
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5))

ggsave(paste0(dir_results, "covariate_balance/weighted_pop_", n_rows, "rows_", modifications, ".png"), weighted_cov_bal_plot)

# print summary statistics for pseudopopulation weights
best_weighted_pseudopop <- merge(ADRD_agg_lagged_trimmed_1_99, subset(gps_for_weighting_list[[best_attempt]],
                                                                      select = c("zip", "year", "capped_stabilized_ipw")),
                                 by = c("zip", "year"))
ess(best_weighted_pseudopop$capped_stabilized_ipw) # 7,889,768

# model outcome from GPS-weighted pseudo-population
formula_expos_only <- as.formula(paste(outcome_name, "~", paste(c("w", strata_vars), collapse = "+", sep = "")))
formula_expos_only_smooth <- as.formula(paste(outcome_name, "~", paste(c("s(w, bs = 'ts')", strata_vars), collapse = "+", sep = "")))

# parametric model (Poisson regression)
bam_exposure_only_capped_weighted <- bam(formula_expos_only,
                                         data = best_weighted_pseudopop,
                                         offset = log(n_persons * n_years),
                                         family = poisson(link = "log"),
                                         weights = capped_stabilized_ipw,
                                         samfrac = 0.05,
                                         chunk.size = 5000,
                                         control = gam.control(trace = TRUE))
summary(bam_exposure_only_capped_weighted)
saveRDS(summary(bam_exposure_only_capped_weighted), file = paste0(dir_results, "parametric_results/bam_capped_weighted_exposure_only_", n_rows, "rows_", modifications, ".rds"))

# semi-parametric model (thin-plate spline)
bam_smooth_exposure_only_capped_weighted <- bam(formula_expos_only_smooth,
                                                data = best_weighted_pseudopop,
                                                offset = log(n_persons * n_years),
                                                family = poisson(link = "log"),
                                                weights = capped_stabilized_ipw,
                                                samfrac = 0.05,
                                                chunk.size = 5000,
                                                control = gam.control(trace = TRUE),
                                                nthreads = n_cores - 1)
png(paste0(dir_results, "semiparametric_results/ERFs/bam_smooth_exposure_only_capped_weighted_", n_rows, "rows_", modifications, ".png"))
plot(bam_smooth_exposure_only_capped_weighted, main = paste0("GPS-Weighted, Capped at 10, Smoothed Poisson regression,\nexposure only (", exposure_name, ")"))
dev.off()
saveRDS(bam_smooth_exposure_only_capped_weighted, file = paste0(dir_results, "semiparametric_results/spline_objects/bam_smooth_exposure_only_capped_weighted_", n_rows, "rows_", modifications, ".rds"))