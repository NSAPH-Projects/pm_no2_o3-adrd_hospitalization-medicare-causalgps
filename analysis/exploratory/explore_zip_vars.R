library(data.table)
library(fst)
library(xtable)
library(ggplot2)
library(reshape2)

# directories for data, code, and results
dir_data <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/data/"
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/code/"
dir_results <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/results/"

# get classifications of variables
source(paste0(dir_code, "analysis/helper_functions.R"))

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

# get distribution of categorical ZIP year variable(s)
for (var in zip_unordered_cat_var_names){
  print(round(prop.table(table(zip_year_data[[var]])), 2))
}
# MW   NE   S    W 
# 0.20 0.19 0.38 0.23

# get mean and SD of quantitative ZIP year covariates
zip_quant_var_distributions <- data.table(Variable = zip_quant_var_names,
                                          Mean = -100, # placeholder value
                                          SD = -100) # placeholder value
for (var in zip_quant_var_names){
  if (var == "summer_tmmx"){
    var_data <- zip_year_data[[var]] - 273.15 # convert from Kelvin to Celsius
  } else var_data <- zip_year_data[[var]]
  zip_quant_var_distributions[Variable == var, `:=`(Mean = round(mean(var_data), 2),
                                                    SD = round(sd(var_data), 2))]
}
print(xtable(zip_quant_var_distributions),
      include.rownames = F)

# get pairwise correlations between ZIP year variables
cor_matrix <- cor(subset(zip_year_data, select = c(zip_expos_names, zip_quant_var_names)))
cor_matrix <- round(cor_matrix, 2)
cor_matrix[upper.tri(cor_matrix)] <- NA
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
ggsave(paste0(dir_results, "exploratory/zip_variables_correlations.png"),
       cor_matrix_plot,
       width = 12,
       units = "in")