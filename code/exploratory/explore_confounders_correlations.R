rm(list = ls())
gc()

##### 1. Setup #####
library(data.table)
library(fst)
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


##### 2. Check correlation between confounders #####

zip_year_data[, `:=`(pm25 = NULL,
                     no2 = NULL,
                     ozone_summer = NULL)]
cor_matrix <- cor(as.matrix(zip_year_data))
cor_matrix <- round(cor_matrix, 2)
cor_matrix[upper.tri(cor_matrix, diag = T)] <- NA # set redundant values to NA; set the diagonal to NA too since it is equal to 1
colnames(cor_matrix) <- c("Summer Humidity",
                          "Summer Temperature",
                          "% Owner Occupied",
                          "Population Density",
                          "Education",
                          "Poverty",
                          "Home Price:Income",
                          "% Black",
                          "% Hispanic",
                          "Smoking Rate",
                          "Body Mass Index")
rownames(cor_matrix) <- colnames(cor_matrix)
cor_matrix <- melt(cor_matrix, na.rm = T)

# generate and save matrix of pairwise correlations
# see http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
cor_matrix_plot <- ggplot(cor_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "red", high = "blue", limit = c(-1,1)) +
  geom_text(aes(Var1, Var2, label = value), color = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank())
cor_matrix_plot
