library(data.table)

# directories for data, code, and results
dir_data <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/data/"
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/code/"
dir_results <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/results/"

# get helpful constants
source(paste0(dir_code, "analysis/constants.R"))

parametric_results_table <- expand.grid(exposure = zip_expos_names, # note: important to NOT call this column "exposure_name", or else won't be able to match with results later
                                        method = c("associational", "weighting", "matching"),
                                        coefficient = -1000.1, # note: important to NOT call this column "coef", or else won't be able to match with results later
                                        se_unadjusted = -1000.1,
                                        se_bootstrap = -1000.1) # these are placeholder numbers (important that they are doubles)

fwrite(parametric_results_table,
       paste0(dir_results, "parametric_results/parametric_results_table.csv"))