### Note: This script should be run after the data have been cleaned, processed, and aggregated data but NOT trimmed (i.e., files #1-5 in the code/aggregation folder have been run) ###

library(data.table)
library(fst)
library(plotrix)
library(ggplot2)

# get directories and classifications of variables
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/git/code/"
source(paste0(dir_code, "constants.R"))


# read in full data
ADRD_agg_lagged <- read_fst(paste0(dir_data, "analysis/ADRD_complete_tv.fst"),
                            as.data.table = TRUE)
setnames(ADRD_agg_lagged,
         old = c("pct_blk", "pct_owner_occ"),
         new = c("prop_blk", "prop_owner_occ"))

# get variables at ZIP year level
zip_year_data <- subset(ADRD_agg_lagged,
                        select = c("zip", "year", zip_expos_names, zip_var_names))
zip_year_data <- unique(zip_year_data, by = c("zip", "year"))

# get summary statistics of ZIP year data
zip_exposure_summary <- expand.grid(Exposure = zip_expos_names,
                                    Mean = -1000.1,
                                    SD = -1000.1,
                                    IQR = -1000.1,
                                    Min = -1000.1,
                                    Fifth = -1000.1,
                                    Median = -1000.1,
                                    NinetyFifth = -1000.1,
                                    Max = -1000.1) # placeholder values; double format
zip_exposure_summary <- as.data.table(zip_exposure_summary)
for (expos in zip_expos_names){
  zip_exposure_summary[Exposure == expos, Mean := mean(zip_year_data[[expos]])]
  zip_exposure_summary[Exposure == expos, SD := sd(zip_year_data[[expos]])]
  zip_exposure_summary[Exposure == expos, IQR := IQR(zip_year_data[[expos]])]
  zip_exposure_summary[Exposure == expos, Min := min(zip_year_data[[expos]])]
  zip_exposure_summary[Exposure == expos, Fifth := quantile(zip_year_data[[expos]], 0.05)]
  zip_exposure_summary[Exposure == expos, Median := median(zip_year_data[[expos]])]
  zip_exposure_summary[Exposure == expos, NinetyFifth := quantile(zip_year_data[[expos]], 0.95)]
  zip_exposure_summary[Exposure == expos, Max := max(zip_year_data[[expos]])]
}
# zip_exposure_summary <- ifelse(zip_exposure_summary > 10, round(zip_exposure_summary, 0),
#                                     round(zip_exposure_summary, 2))
fwrite(zip_exposure_summary, paste0(dir_results, "exploratory/zip_exposure_summary.csv"))

# get pairwise correlations of ZIP year exposures
# > cor(zip_year_data$pm25, zip_year_data$no2)
# [1] 0.3991182
# > cor(zip_year_data$pm25, zip_year_data$ozone_summer)
# [1] 0.2499757
# > cor(zip_year_data$no2, zip_year_data$ozone_summer)
# [1] 0.2503145
