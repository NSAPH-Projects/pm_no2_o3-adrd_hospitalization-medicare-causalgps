# get directories and classifications of variables
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/git/code/"
source(paste0(dir_code, "constants.R"))

# initiate table for regression results (coefficient for exposure)
cat(paste("Exposure", "Method", "Coefficient", sep = ","),
    sep = "\n",
    file = paste0(dir_results, "parametric_results/coef_for_exposure.txt"),
    append = TRUE)

# initiate table for m-out-of-n bootstrap standard errors (of coefficient for exposure)
cat(paste("Exposure", "Method", "BootSE", sep = ","),
    sep = "\n",
    file = paste0(dir_results, "parametric_results/boot_SE.txt"),
    append = TRUE)
