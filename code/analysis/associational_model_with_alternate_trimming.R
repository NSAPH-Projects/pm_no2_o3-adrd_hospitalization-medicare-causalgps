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


### Sensitivity analysis: thin-plate spline outcome model ###

# function to fit model
bam_smooth_trimmed <- function(exposure_name, dataset, trim_level){
  cl <- parallel::makeCluster(n_cores, type = "PSOCK")
  bam_smooth_associational <- 
    bam(as.formula(paste("Y ~", paste(c("s(w, bs = 'cr', k = 4, fx = TRUE)",
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
  
  # estimate counterfactual for every year-zip-strata, calculate ATE
  potential_data <- copy(dataset)
  data_prediction <- 
    rbindlist(lapply(seq(min(dataset$w), 
                         max(dataset$w), 
                         length.out = 100), function(pot_exp) {
                           
                           # Get potential data if all had same potential exposure
                           potential_data[, w := pot_exp]
                           return(data.table(name = exposure_name,
                                             w = pot_exp,
                                             ate = mean(predict(bam_smooth_associational, 
                                                                newdata = potential_data, 
                                                                type = "response"))))
                         }))
  plot(I(1e5*ate)~w,data_prediction, type = 'l')
  save(data_prediction, file = paste0(dir_results, "semiparametric_results/", exposure_name, "_trim_", trim_level[1], "_", trim_level[2], "_assoc_smooth.rda"))
}

# to do: run the spline for trim_level_1 and trim_level_2 over zip_expos_names?
