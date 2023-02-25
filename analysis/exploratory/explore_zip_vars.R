library(ggplot2)
library(reshape2)

zip_year_data <- subset(ADRD_agg_lagged,
                        select = c("zip", "year", zip_expos_names, zip_var_names))
zip_year_data <- unique(zip_year_data, by = c("zip", "year"))

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
# > cor(zip_year_data$ozone_summer, zip_year_data$pm25)
# [1] 0.2499757
# > cor(zip_year_data$no2, zip_year_data$pm25)
# [1] 0.3991182
# > cor(zip_year_data$no2, zip_year_data$ozone_summer)
# [1] 0.2503145

temp <- cor(subset(zip_year_data, select = c(zip_expos_names, zip_quant_var_names)))
temp <- round(temp, 2)
temp[upper.tri(temp)] <- NA
temp <- melt(temp, na.rm = T)

p <- ggplot(temp, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "red", high = "blue", limit = c(-1,1)) +
  geom_text(aes(Var1, Var2, label = value), color = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank())
ggsave(paste0(dir_results, "exploratory/zip_variables_correlations.png"), p, width = 12, units = "in")