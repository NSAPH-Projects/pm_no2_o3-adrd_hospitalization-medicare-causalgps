### Note: This script should be run after the data have been cleaned, processed, and aggregated data but NOT trimmed (i.e., files #1-5 in the code/aggregation folder have been run) ###

library(data.table)
library(fst)

# get directories and classifications of variables
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/git/code/"
source(paste0(dir_code, "constants.R"))


# read in full data (observations are aggregated patient units)
ADRD_agg_lagged <- read_fst(paste0(dir_data, "analysis/ADRD_complete_tv.fst"),
                            as.data.table = TRUE)
setnames(ADRD_agg_lagged,
         old = c("pct_blk", "pct_owner_occ"),
         new = c("prop_blk", "prop_owner_occ"))

# calculate pairwise correlations of exposures at aggregated patient unit level
pairwise_correlations <- data.table(Exposure1 = zip_expos_names,
                                    Exposure2 = c(zip_expos_names[2:length(zip_expos_names)], zip_expos_names[1]),
                                    Correlation = NA)
for (i in 1:nrow(pairwise_correlations)){
  pairwise_correlations$Correlation[i] <- cor(ADRD_agg_lagged[[pairwise_correlations$Exposure1[i]]],
                                              ADRD_agg_lagged[[pairwise_correlations$Exposure2[i]]])
}
pairwise_correlations
