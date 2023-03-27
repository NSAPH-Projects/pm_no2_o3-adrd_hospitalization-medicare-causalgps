##### Setup #####

library(data.table)
library(fst)
library(mgcv)

# directories for data, code, and results
dir_data <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/data/"
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/code/"
dir_results <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/results/"

# set exposure
exposure_name <- "pm25"

# parameters for this computing job
n_cores <- 1 # 48 is max of fasse partition, 64 is max of fasse_bigmem partition

# get m and n values for this m-out-of-n bootstrap
if (exposure_name == "pm25"){
  m_boot <- 2964
} else if (exposure_name == "no2"){
  m_boot <- 2990
} else if (exposure_name == "ozone_summer"){
  m_boot <- 2937
} else message("'exposure_name' must be 'pm25', 'no2', or 'ozone_summer'")

# get helpful constants
source(paste0(dir_code, "analysis/constants.R"))

# get m out of n bootstrap sample
boot_sample_number <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

boot_zips <- read.csv(paste0(dir_data, "analysis/",
                             exposure_name, "/",
                             "boot_zips/",
                             m_boot, "zips_1000replicates.csv"))
boot_zips <- boot_zips[, boot_sample_number]

zip_year_data_with_strata <- read_fst(paste0(dir_data, "analysis/", exposure_name, "/zip_year_data_with_strata_trimmed_0.05_0.95.fst"), as.data.table = T)
boot_zip_year_data_with_strata <- zip_year_data_with_strata[zip %in% boot_zips]
rm(zip_year_data, zip_year_data_with_strata)

# make sure categorical variables are factors
boot_zip_year_data_with_strata[, `:=`(zip = as.factor(zip),
                                      year = as.factor(year),
                                      age_grp = as.factor(age_grp),
                                      sex = as.factor(sex),
                                      race = as.factor(race),
                                      dual = as.factor(dual))]

##### Run model (Poisson regression) #####

# run parametric outcome model
bam_associational <- bam(as.formula(paste("Y ~", paste(c("w",
                                                         strata_vars,
                                                         zip_var_names),
                                                       collapse = "+", sep = ""))),
                         data = boot_zip_year_data_with_strata,
                         offset = log(n_persons * n_years),
                         family = poisson(link = "log"),
                         samfrac = 0.05,
                         chunk.size = 5000,
                         control = gam.control(trace = TRUE),
                         nthreads = n_cores)

# result to be bootstrapped: coefficient for exposure
coef_for_exposure <- summary(bam_associational)$p.table["w", "Estimate"]
result <- data.table(boot_sample_number = boot_sample_number,
                     coef_for_exposure = coef_for_exposure)
fwrite(result,
       paste0(dir_results, "bootstrap/",
              exposure_name, "/",
              "associational/",
              m_boot, "zips/",
              "replicate_", boot_sample_number, ".csv"))
