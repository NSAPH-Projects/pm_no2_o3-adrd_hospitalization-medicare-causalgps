library(mgcv)
library(ggplot2)
# > names(ADRD_agg)
# [1] "AD_zip"           "ADRD_year"          "ffs_entry_year"    
# [4] "ADRD_age"           "entry_sex"          "entry_race"        
# [7] "exit_dual"            "n_persons"          "n_ADRD_hosp"        
# [10] "n_years"            "conf_year"          "mean_bmi"          
# [13] "smoke_rate"         "hispanic"           "pct_blk"           
# [16] "medhouseholdincome" "medianhousevalue"   "poverty"           
# [19] "education"          "popdensity"         "pct_owner_occ"     
# [22] "summer_tmmx"        "winter_tmmx"        "summer_rmax"       
# [25] "winter_rmax"        "city"               "statecode"         
# [28] "latitude"           "longitude"          "pm25"              
# [31] "no2"                "ozone_summer"       "ozone_winter"      
# [34] "ozone_fall"         "ozone_spring"       "tmmx"              
# [37] "rmax"               "pr"  

# Covariate balance
bal <- cor(ADRD_agg[, pm25], 
           ADRD_agg[, .(ADRD_year, ADRD_age, sexM, any_dual,
                        white=race_cat=="white", black=race_cat=="black", other=race_cat=="other", 
                        hisp=race_cat=="hisp", n_am_nat=race_cat=="n_amer_native", asian=race_cat=="asian",
                        mean_bmi, smoke_rate, hispanic, pct_blk, PIR,
                        education, popdensity, pct_owner_occ)])
baldf <- data.frame(n = colnames(bal)[order(abs(bal))], 
                    y = 1:length(bal),
                    x = sort(abs(bal)))
baldf$n <- factor(baldf$n, levels = colnames(bal)[order(abs(bal))],
                  labels = colnames(bal)[order(abs(bal))])
ggplot(baldf) +
  geom_point(aes(x = x, y = y)) +
  geom_line(aes(x = x, y = y)) +
  geom_vline(xintercept = 0.1, linetype = 2) +
  theme_bw() +
  scale_y_continuous(breaks = 1:length(bal), labels = baldf$n)


# Calculate GPS 
gps <- bam(pm25 ~ factor(ADRD_year) + factor(any_dual) +
             ADRD_age + factor(sexM) + factor(race_cat) +
             mean_bmi + smoke_rate + hispanic + pct_blk +
             PIR + poverty +
             education + popdensity + pct_owner_occ,
           data = ADRD_agg)
summary(gps)
# Stabilized IPWs
ADRD_agg[, ipw_pm25 := dnorm(pm25, mean(pm25), sd(pm25)) /
           dnorm(pm25, gps$fitted.values, sqrt(mean(gps$residuals^2)))]
ADRD_agg[ipw_pm25 > 10, ipw_pm25 := 10]
# Weighted balance
bal2 <- cov.wt(ADRD_agg[, .(pm25, 
                            ADRD_year, ADRD_age, sexM, any_dual,
                            white=race_cat=="white", black=race_cat=="black", 
                            other=race_cat=="other", hisp=race_cat=="hisp", 
                            n_am_nat=race_cat=="n_amer_native", asian=race_cat=="asian",
                            mean_bmi, smoke_rate, hispanic, pct_blk, PIR,
                            education, popdensity, pct_owner_occ)],
               ADRD_agg[, ipw_pm25], cor = TRUE)$cor[,1]
bal2df <- merge(baldf, data.frame(n = names(bal2), x2 = abs(bal2)), by = c("n"))
ggplot(bal2df) +
  geom_point(aes(x = x, y = y)) +
  geom_line(aes(x = x, y = y)) +
  geom_point(aes(x = x2, y = y), color = 2) +
  geom_vline(xintercept = 0.1, linetype = 2) +
  theme_bw() +
  scale_y_continuous(breaks = 1:length(bal), labels = baldf$n)

# Model without weights
m1 <- bam(n_ADRDhosp ~ offset(log(n_persons * n_years)) +
            pm25 + no2 + ozone_summer + tmmx + rmax +
            factor(ffs_entry_year) + factor(n_years) + factor(any_dual) +
            ADRD_age + factor(sexM) + factor(race_cat) +
            mean_bmi + smoke_rate + hispanic + pct_blk +
            PIR + poverty +
            education + popdensity + pct_owner_occ,
          data = ADRD_agg,
          family = poisson, 
          samfrac = 0.05, chunk.size = 5000,
          control = gam.control(trace = TRUE))
summary(m1)

# Causal inference model
m2 <- bam(n_ADRDhosp ~ offset(log(n_persons * n_years)) +
            pm25 + no2 + ozone_summer + tmmx + rmax +
            factor(ffs_entry_year) + factor(n_years) + factor(any_dual) +
            ADRD_age + factor(sexM) + factor(race_cat) +
            mean_bmi + smoke_rate + hispanic + pct_blk +
            PIR + poverty +
            education + popdensity + pct_owner_occ,
          data = ADRD_agg,
          weights = ADRD_agg[, ipw_pm25],
          family = poisson, 
          samfrac = 0.05, chunk.size = 5000,
          control = gam.control(trace = TRUE))
summary(m2)

# HR of IQR increase in exposure
exp(m1$coefficients[2]*IQR(ADRD_agg$pm25))

exp(m2$coefficients[2]*IQR(ADRD_agg$pm25))
