library(data.table)
library(fst)
library(plotrix)
library(ggplot2)

# directories for data, code, and results
dir_data <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/data/"
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/code/"
dir_results <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/results/"

# get classifications of variables
source(paste0(dir_code, "analysis/constants.R"))

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

# get mean, SD, IQR of ZIP year exposures
# > mean(zip_year_data$pm25)
# [1] 9.951262
# > sd(zip_year_data$pm25)
# [1] 3.223463
# > mean(zip_year_data$no2)
# [1] 16.91286
# > sd(zip_year_data$no2)
# [1] 9.430883
# > mean(zip_year_data$ozone_summer)
# [1] 45.99236
# > sd(zip_year_data$ozone_summer)
# [1] 7.514933
# > IQR(zip_year_data$pm25)
# [1] 4.159233
# > IQR(zip_year_data$no2)
# [1] 12.01235
# > IQR(zip_year_data$ozone_summer)
# [1] 9.800937

# get pairwise correlations of ZIP year exposures
# > cor(zip_year_data$pm25, zip_year_data$no2)
# [1] 0.3991182
# > cor(zip_year_data$pm25, zip_year_data$ozone_summer)
# [1] 0.2499757
# > cor(zip_year_data$no2, zip_year_data$ozone_summer)
# [1] 0.2503145
