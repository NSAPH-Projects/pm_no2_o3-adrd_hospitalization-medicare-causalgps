# Summary of the run:
#
# Objective: Generating pseudo population based on matching approach for the 
# following parameters.
#

#rm(list = ls())
#gc()

setwd("/n/dominici_nsaph_l3/Lab/projects/nkhoshnevis_pm_no2_o3-adrd_hosp-medicare-causalgps/pm_no2_o3-adrd_hospitalization-medicare-causalgps/analysis")

trim_mat <- matrix(c(0.01, 0.05, 0.1, 0.15,  0.99, 0.95, 0.9, 0.85),
                   ncol = 2, byrow = FALSE)
delta_n_vector <- c(0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0)
all_runs <- expand.grid((1:nrow(trim_mat)), delta_n_vector)
run_name <- paste0("1221_with_noise_", 1:nrow(all_runs))

args <-  commandArgs(trailingOnly = TRUE)
job_array_number <- as.integer(args[1])
print(paste0("This is from task ID: ", job_array_number))

set.seed(4987)
seed_vector <- sample(100000:1000000, size = nrow(all_runs))
seed_val <- seed_vector[job_array_number]

# parameters

run_number <- run_name[job_array_number]
exposure_name <- "pm25"
outcome_name <- "n_hosp"
n_cores <- 36
trim_quantiles <- c(trim_mat[all_runs[1][[1]][job_array_number], 1], 
                    trim_mat[all_runs[1][[1]][job_array_number], 2])
delta_n <- all_runs[2][[1]][job_array_number]


output_folder <- paste0("cgps", 
                        "_exp_", exposure_name,
                        "_out_", outcome_name, 
                        "_trimmed_",as.integer(trim_quantiles[1]*100),"_",
                        as.integer(trim_quantiles[2]*100),
                        "_delta_", delta_n, 
                        "_run_", run_number)

if (!dir.exists(output_folder)) {dir.create(output_folder)}

print(output_folder)

# 0. Setup ---------------------------------------------------------------------
# Un-comment the folloowing line to install the CausalGPS package

# options(timeout=1000000)
# devtools::install_github("NSAPH-Software/CausalGPS",
#                          ref="develop")

library(data.table)
library(fst)
library(CausalGPS)
library(mgcv)
library(ggplot2)
library(tidyr)

# directory
readRenviron(".Renviron")
dir_proj <- Sys.getenv("PROJ_PATH_pm_no2_o3_adrd_hosp_medicare_causalgps")

# Load helper functions
source(paste0(dir_proj, "code/analysis/helper_functions.R"))
# ------------------------------------------------------------------------------

# 1. Read and preprocess data --------------------------------------------------

ADRD_agg_lagged <- fst::read_fst(paste0(dir_proj,
                                        "data/analysis/ADRD_complete_tv.fst"),
                                 as.data.table = TRUE)
data.table::setnames(ADRD_agg_lagged, old = c("pct_blk", "pct_owner_occ"),
                                      new = c("prop_blk", "prop_owner_occ"))

ADRD_agg_lagged[, `:=`(zip = as.factor(zip), 
                       year = as.factor(year), 
                       cohort = as.factor(cohort), 
                       age_grp = as.factor(age_grp), 
                       sex = as.factor(sex), 
                       race = as.factor(race), 
                       dual = as.factor(dual))]

# Get data for exposure, outcome, and covariates of interest -----

#exposure_name <- "pm25"
other_expos_names <- zip_expos_names[zip_expos_names != exposure_name]
#outcome_name <- "n_hosp"
ADRD_agg_lagged_subset <- subset(ADRD_agg_lagged,
                                 select = c(exposure_name, 
                                            outcome_name, 
                                            other_expos_names, 
                                            zip_var_names, 
                                            indiv_var_names, 
                                            offset_var_names))


# add some noise to the exposure value ---------------------------
max_noise <- 0.025
pm_noise <- runif(nrow(ADRD_agg_lagged_subset))*max_noise - 
  runif(nrow(ADRD_agg_lagged_subset))*max_noise
ADRD_agg_lagged_subset$pm_noise <- pm_noise
ADRD_agg_lagged_subset[, `:=`(pm25 = pm25*(1+pm_noise))]

for (var in c(zip_unordered_cat_var_names, indiv_unordered_cat_var_names)){
  ADRD_agg_lagged_subset[[var]] <- as.factor(ADRD_agg_lagged_subset[[var]])
}

# Trim exposures outside the 1st and 99th percentiles, for all analyses 
# (associational and causal) ----

trimmed_data <- quantile(ADRD_agg_lagged_subset[[exposure_name]], trim_quantiles)
rows_within_range <- ADRD_agg_lagged_subset[[exposure_name]] >= trimmed_data[1] & 
                     ADRD_agg_lagged_subset[[exposure_name]] <= trimmed_data[2]
ADRD_agg_lagged_trimmed <- ADRD_agg_lagged_subset[rows_within_range, ]

print(paste0(" Number of original data: ", nrow(ADRD_agg_lagged_subset),
             " Number of data after trimming: ", 
             nrow(ADRD_agg_lagged_trimmed)))

n_rows <- nrow(ADRD_agg_lagged_trimmed)

# ------------------------------------------------------------------------------

# 2. Prepare data to the analyses ----------------------------------------------

# Option 1: all data
selected_rows <- 1:n_rows 

# Option 2: subset of data
#set.seed(3986)
#selected_rows <- sample(seq(1,n_rows), 100000, replace = FALSE, prob = NULL)


# TODO: it should also work with data.table. If not, we need to add it. 
Y <- as.data.frame(ADRD_agg_lagged_trimmed)[selected_rows, outcome_name]
w <- as.data.frame(ADRD_agg_lagged_trimmed)[selected_rows, exposure_name]
c_data <- as.data.frame(subset(ADRD_agg_lagged_trimmed[selected_rows,],
                        select = c(other_expos_names, zip_var_names)))

# ------------------------------------------------------------------------------

# 3. CausalGPS runs ------------------------------------------------------------

# parameters for this computing job
# n_cores <- 48 # 48

# to be used in names of output files, to record how you're tuning the models
# modifications <- "bin_age_tv_trimmed" 

set_logger(logger_file_path = paste0(output_folder,"/CausalGPS.log"),
           logger_level = "TRACE")

set.seed(seed_val)
matched_pop_subset <- generate_pseudo_pop(Y,
                                          w,
                                          c_data,
                                          ci_appr = "matching",
                                          gps_model = "parametric",
                                          use_cov_transform = TRUE,
                                          transformers = list("pow2",
                                                              "pow3",
                                                              "sqrt", 
                                                              "log_nonneg", 
                                                              "logit_nonneg"),
                                          sl_lib = c("m_xgboost"),
                                          #bin_seq = seq(4,16,0.2),
                                          params = list(xgb_nrounds = seq(10, 60),
                                                        xgb_eta = seq(0.1, 0.4, 0.01)),
                                          nthread = n_cores,
                                          covar_bl_method = "absolute",
                                          covar_bl_trs = 0.1,
                                          covar_bl_trs_type = "maximal",
                                          optimized_compile = TRUE,
                                          trim_quantiles = c(0,1),
                                          max_attempt = 40,
                                          matching_fun = "matching_l1",
                                          delta_n = delta_n,
                                          scale = 1)

pdf(paste0(output_folder,"/PseudoPop.pdf"), width = 8, height = 8)
plot(matched_pop_subset)
dev.off()

saveRDS(matched_pop_subset,
        file = paste0(output_folder,"/matched_pop_subset.rds"))

writeLines(capture.output(sessionInfo()),  
           paste0(output_folder,"/sessioninfo.txt"))

writeLines(capture.output(summary(matched_pop_subset)),
           paste0(output_folder, "/summary.txt"))

# ------------------------------------------------------------------------------
