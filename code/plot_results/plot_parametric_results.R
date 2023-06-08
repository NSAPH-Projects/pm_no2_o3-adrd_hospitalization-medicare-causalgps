library(data.table)
library(ggplot2)

# get directories and classifications of variables
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/git/code/"
source(paste0(dir_code, "constants.R"))

# get hazard ratio results from all models
parametric_results_table <- fread(paste0(dir_results, "parametric_results/all_HR_results.csv"))

# plot hazard ratios with 95% confidence intervals
ggplot(parametric_results_table, aes(x = Exposure,
                                          y = HR,
                                          color = Method,
                                          shape = Method)) +
  geom_hline(yintercept = 1, linewidth = 1, linetype = 2) +
  geom_point(position=position_dodge(0.5), size = 4) +
  geom_errorbar(aes(ymin = HR_95CI_lower,
                    ymax = HR_95CI_upper),
                position=position_dodge(0.5), width = 0.4,
                linewidth = 1) +
  ylab("Hazard Ratio per IQR Increase in Exposure") +
  theme_minimal(base_size = 24) +
  theme(legend.key.width = unit(1, "cm"), 
        legend.key.height = unit(1, "cm"),
        legend.box.spacing = unit(-1, "cm"),
        legend.position = "bottom",
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(color = "", shape = "", x = "") +
  scale_x_discrete(labels = parse(text = c("PM[2.5]", "NO[2]", "Summer~Ozone")))

# # save plot
# ggsave(paste0(dir_results, "parametric_results/all_parametric_results.png"),
#        p,
#        width = 569,
#        height = 368,
#        units = "px")
