library(ggplot2)
library(data.table)

# get directories
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/git/code/"
source(paste0(dir_code, "constants.R"))

res <-
  rbindlist(lapply(c("pm25", "no2", "ozone_summer"), function(e) {
    rbindlist(lapply(c("gpsmatching", "gpsweighting", "assoc"), function(m) {
      if (file.exists(paste0(dir_results, e, "_", m, "_smooth.rda"))) {
        load(paste0(dir_results, e, "_", m, "_smooth.rda"))
        setDT(data_prediction)
        data_prediction[, model := m]
      }
    }))
  }))


res$model <- factor(res$model, c("assoc", "gpsweighting", "gpsmatching"),
                    c("Poisson Regression", "GPS Weighting", "GPS Matching"))
res$name <- factor(res$name, c("pm25", "no2", "ozone_summer"),
                   c("PM[2.5]", "NO[2]", "Summer~Ozone"))

ggplot(res, aes(x = w, y = ate * 1e5, color = model, linetype = model)) +
  geom_line(linewidth = 2) +
  facet_grid(~name, scales = "free_x", labeller = label_parsed) +
  theme_minimal(base_size = 24) +
  guides(color = guide_legend(nrow = 3)) +
  theme(legend.key.width = unit(2, "cm"), 
        legend.key.height = unit(1, "cm"),
        legend.position = c(.88, .2),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"), 
        strip.placement = "outside") +
  labs(x = "Annual average exposures", y = "Hospitalizations with ADRD\nper 100,000 beneficiaries",
       color = "", linetype = "")
