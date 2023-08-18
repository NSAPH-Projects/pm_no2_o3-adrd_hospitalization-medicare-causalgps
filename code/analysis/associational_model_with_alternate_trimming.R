##### Note: For NO2, the bam model may require 60 GB to run #####

rm(list = ls())
gc()

##### 1. Setup #####
library(data.table)
library(fst)
library(mgcv)

# get directories and classifications of variables
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/git/code/"
source(paste0(dir_code, "constants.R"))

# set parameters for this computing job
n_cores <- 4
trim_level_1 <- c(0.03, 0.97)
trim_level_2 <- c(0.01, 0.99)

# read in full data
ADRD_agg_lagged <- read_fst(paste0(dir_data, "analysis/ADRD_complete_tv.fst"),
                            as.data.table = TRUE)
setnames(ADRD_agg_lagged,
         old = c("pct_blk", "pct_owner_occ"),
         new = c("prop_blk", "prop_owner_occ"))
ADRD_agg_lagged[, sex := ifelse(sex == 2, 0, # in raw data, female is coded as 2; recode as 0
                                ifelse(sex == 1, 1, NA))] # in raw data, male is coded as 1; keep this; there should be no NAs
ADRD_agg_lagged_subset <- subset(ADRD_agg_lagged, select = c(zip_expos_names,
                                                             outcome_name,
                                                             zip_var_names,
                                                             "zip",
                                                             indiv_var_names,
                                                             offset_var_names))
for (var in c("zip",
              zip_unordered_cat_var_names,
              indiv_unordered_cat_var_names)){
  ADRD_agg_lagged_subset[[var]] <- as.factor(ADRD_agg_lagged_subset[[var]])
}

# get dataset of 1 observation per ZIP, year; use this to trim exposures
zip_year_data <- subset(ADRD_agg_lagged_subset, select = c("zip",
                                                           "year",
                                                           zip_expos_names))
zip_year_data <- unique(zip_year_data, by = c("zip", "year"))

# function to trim dataset for any exposure and trim level
trim_exposure <- function(exposure_name, trim_level){
  trim <- quantile(zip_year_data[[exposure_name]], trim_level)
  rows_within_range <- ADRD_agg_lagged_subset[[exposure_name]] >= trim[1] & ADRD_agg_lagged_subset[[exposure_name]] <= trim[2]
  ADRD_agg_lagged_subset_trimmed <- ADRD_agg_lagged_subset[zip_year_rows_within_range, ]
  colnames(ADRD_agg_lagged_subset_trimmed)[colnames(ADRD_agg_lagged_subset_trimmed) == exposure_name] <- "w"
  colnames(ADRD_agg_lagged_subset_trimmed)[colnames(ADRD_agg_lagged_subset_trimmed) == outcome_name] <- "Y"
  return(ADRD_agg_lagged_subset_trimmed)
}

# function to run model (Poisson regression)
bam_trimmed <- function(exposure_name, dataset, trim_level){
  cl <- parallel::makeCluster(n_cores, type = "PSOCK")
  bam_associational <- bam(as.formula(paste("Y ~", paste(c("w",
                                                           strata_vars,
                                                           zip_var_names),
                                                         collapse = "+", sep = ""))),
                           data = dataset,
                           offset = log(n_persons * n_years),
                           family = poisson(link = "log"),
                           samfrac = 0.05,
                           chunk.size = 5000,
                           control = gam.control(trace = TRUE),
                           nthreads = n_cores,
                           cluster = cl)
  parallel::stopCluster(cl)
  cat(paste(exposure_name, "Poisson Regression", paste(c("trim", trim_level)), bam_associational$coefficients["w"], sep = ","),
      sep = "\n",
      file = paste0(dir_results, "parametric_results/coef_for_exposure_with_alternate_trimming.txt"),
      append = TRUE)
}

# run the analyses
for (exposure_name in zip_expos_names){
  bam_trimmed(exposure_name, dataset, trim_level_1)
}

for (exposure_name in zip_expos_names){
  bam_trimmed(exposure_name, dataset, trim_level_2)
}

# convert regression coefficient to hazard ratio
coef <- fread(paste0(dir_results, "parametric_results/coef_for_exposure_with_alternate_trimming.txt"))
zip_exposure_summary <- fread(paste0(dir_results, "exploratory/zip_exposure_summary.csv"))
iqr <- zip_exposure_summary[, .(Exposure, IQR)] # 4.159 micrograms/m^2 for PM2.5, 12.012 ppb for NO2, 9.801 ppb for summer ozone
coef <- merge(coef, iqr, by = "Exposure")
coef[, HR := exp(IQR * Coefficient)]
coef
