##### Setup #####

# devtools::install_github("fasrc/CausalGPS", ref="develop")
library(data.table)
library(fst)
library(mgcv)

# directories for data, code, and results
dir_data <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/data/"
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/code/"
dir_results <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/results/"

# set exposure
exposure_name <- "pm25"

# set parameters for this computing job
n_cores <- 16

# get data and helpful functions
source(paste0(dir_code, "analysis/helper_functions.R"))
zip_year_data_with_strata <- read_fst(paste0(dir_data, "analysis/",
                                             exposure_name, "/",
                                             "zip_year_data_with_strata_trimmed_0.05_0.95.fst"),
                                      as.data.table = T)

# make sure categorical variables are factors
zip_year_data_with_strata[, `:=`(zip = as.factor(zip),
                                 year = as.factor(year),
                                 age_grp = as.factor(age_grp),
                                 sex = as.factor(sex),
                                 race = as.factor(race),
                                 dual = as.factor(dual))]

# run model (Poisson regression)
cl <- parallel::makeCluster(n_cores, type = "PSOCK")
bam_associational <- bam(as.formula(paste("Y ~", paste(c("w",
                                                         strata_vars,
                                                         zip_var_names),
                                                       collapse = "+", sep = ""))),
                         data = zip_year_data_with_strata,
                         offset = log(n_persons * n_years),
                         family = poisson(link = "log"),
                         samfrac = 0.05,
                         chunk.size = 5000,
                         control = gam.control(trace = TRUE),
                         nthreads = n_cores,
                         cluster = cl)
parallel::stopCluster(cl)

# save results
saveRDS(summary(bam_associational),
        file = paste0(dir_results, "parametric_results/",
                      exposure_name, "/",
                      "associational/",
                      "bam_parametric_associational_", nrow(zip_year_data_with_strata), "rows.rds"))

# sensitivity analysis: thin-plate spline model
cl <- parallel::makeCluster(n_cores, type = "PSOCK")
bam_smooth_associational <- bam(as.formula(paste("Y ~", paste(c("s(w, bs = 'ts')",
                                                         strata_vars,
                                                         zip_var_names),
                                                       collapse = "+", sep = ""))),
                         data = zip_year_data_with_strata,
                         offset = log(n_persons * n_years),
                         family = poisson(link = "log"),
                         samfrac = 0.05,
                         chunk.size = 5000,
                         control = gam.control(trace = TRUE),
                         nthreads = n_cores,
                         cluster = cl)
parallel::stopCluster(cl)

# save results for sensitivity analysis
png(paste0(dir_results, "semiparametric_results/ERFs/",
           exposure_name, "/",
           "associational/",
           "bam_smooth_associational_", nrow(zip_year_data_with_strata), "rows.png"))
plot(bam_smooth_associational,
     main = paste0("Associational, Smoothed Poisson regression\nExposure: ", exposure_name))
dev.off()