###############################################################################
# Project: Causal effect of air pollution on first ADRD hospitalization       #
# Code: non-causal linear and smooth effects of PM2.5, NOx                    #
# Input: "ADRD_aggregated_entry_to_followup_merged.fst"                       #
# Author: Daniel Mork                                                         #
# Date: 2021-12-28                                                            #
###############################################################################


######################## ? Questions/considerations ############################
#' 1. Age bins produce very low counts of ADRD events, need for zero-inflated
#'    type model?
#'    Ans:
#' 2. Very different results from Liu/Wu neuro paper?
#'    Ans:
#' 3. Control for entry year and followup year?
#'    Ans:
#' 4. Connect zipcode census data to entry or followup year?
#'    Ans:
#' 
############################# Things to try ####################################
#' Try: larger age bins (2 yr, 5 yr) to reduce zeros in ADRD event count
#' Outcome:
#' Try: categorical age effects (vs. 2nd degree poly)
#' Outcome:
#' Try: cox model
#' Outcome:


rm(list = ls())
gc()
############################# 0. Setup ########################################
library(NSAPHutils)
library(data.table)
library(fst)
library(mgcv)
library(parallel)
setDTthreads(threads = 16)

dir_data <- "/nfs/home/D/dam9096/shared_space/ci3_analysis/dmork/Data/Causal_ADRD/"

############################# 1. Load data #####################################
ADRDdat <- read_fst(paste0(dir_data, "ADRD_aggregated_entry_to_followup_merged.fst"), as.data.table = TRUE)
names(ADRDdat)
# [1] "followupZip"  "followupYear" "sex"          "race"         "dual"        
# [6] "entryYear"    "age"          "nADRDEvents"  "nDeaths"      "N"           
# [11] "nYears"       "followupAge"  "avg_age"      "med_age"      "bmi"         
# [16] "smoke"        "hispanic"     "pct_blk"      "income"       "value"       
# [21] "poverty"      "education"    "popdensity"   "ownerocc"     "summer_tmmx" 
# [26] "summer_rmax"  "winter_tmmx"  "winter_rmax"  "death_prob"   "pm25"        
# [31] "no2"          "ozone"        "nox"          "city"         "state"       
# [36] "latitude"     "longitude"    "region"      
ADRDdat <- ADRDdat[complete.cases(ADRDdat)]
ADRDdat[, .(sum(N), sum(nYears * N), sum(nADRDEvents))]
# 71937782 unique individuals
# 599875618 person-years
# 8074368 ADRD events
ADRDdat[, .(sum(N),  sum(nADRDEvents)), by = age][order(age)]
# age       V1     V2
# 1:  65 40114996 892550
# 2:  66  2080559 150803
# 3:  67  1978438 170591
# 4:  68  1900711 185963
# 5:  69  1896242 212721
# 6:  70  1774386 226166
# 7:  71  1755818 249338
# 8:  72  1727867 270874
# 9:  73  1650972 283644
# 10:  74  1604433 300031
# 11:  75  1568631 319621
# 12:  76  1472174 319310
# 13:  77  1395870 322341
# 14:  78  1363764 332659
# 15:  79  1244236 319814
# 16:  80  1073413 289513
# 17:  81  1021621 285541
# 18:  82   913860 264524
# 19:  83   823482 245585
# 20:  84   748430 228682
# 21:  85   676465 209986
# 22:  86   583223 182923
# 23:  87   512084 161990
# 24:  88   425173 133196
# 25:  89   362288 112599
# 26:  90   300622  92937
# 27:  91   249289  75174
# 28:  92   200547  59477
# 29:  93   154723  44310
# 30:  94   120708  33257
# 31:  95    90178  24031
# 32:  96    67656  16948
# 33:  97    51724  11853
# 34:  98   165425  18920
# 35:  99       87      9
# 36: 100       68      4
# 37: 101       40      2
# 38: 102       30      2
# 39: 103       17      0
# 40: 104       19      1
# 41: 105       15      0
# 42: 106       19      1
# 43: 107        8      0
# 44: 108       10      0
# 45: 109        4      0
# 46: 110        5      0
# 47: 111        4      0
# 48: 112        5      0
# 49: 113        3      0
# 50: 114       43      0
# age       V1     V2
ADRDdat[, .(sum(N),  sum(nADRDEvents)), by = sex]
# sex       V1      V2
# 1:   1 31988374 2513533
# 2:   2 40082011 4534358
ADRDdat[, .(sum(N),  sum(nADRDEvents)), by = dual]
# dual       V1      V2
# 1:    1  8309902 1241636
# 2:    0 63760483 5806255
ADRDdat[, .(sum(N),  sum(nADRDEvents)), by = race]
# race       V1      V2
# 1:  blk  6529151  758486
# 2:  his  1459977  127921
# 3:  wht 61200871 6029544
# 4:  oth  2880386  131940
ADRDdat[nADRDEvents == 0, .N]
# [1] 14,903,275 out of 21+M zero events. need zero-inflated model? smaller age bins?


########################## ** Poisson models ** ################################
########################## 2. Linear effects ###################################
m1 <- glm(nADRDEvents ~ 
            offset(log(I(nYears * N))) + # offset total person-time at-risk
            pm25 + 
            nox + 
            poly(followupAge, 2) +
            as.factor(region) + 
            as.factor(followupYear) +
            as.factor(sex) + as.factor(race) + as.factor(dual) +
            hispanic + pct_blk +
            bmi + smoke + income + value + poverty +
            education + popdensity + ownerocc +
            summer_tmmx + summer_rmax + winter_tmmx + winter_rmax,
          data = ADRDdat, model = F,
          family = poisson(link = "log")) # print diagnostic output
summary(m1)
# Call:
#   glm(formula = nADRDEvents ~ offset(log(I(nYears * N))) + pm25 + 
#         nox + poly(followupAge, 2) + as.factor(region) + as.factor(followupYear) + 
#         as.factor(sex) + as.factor(race) + as.factor(dual) + hispanic + 
#         pct_blk + bmi + smoke + income + value + poverty + education + 
#         popdensity + ownerocc + summer_tmmx + summer_rmax + winter_tmmx + 
#         winter_rmax, family = poisson(link = "log"), data = ADRDdat, 
#       model = F)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -4.8855  -0.7176  -0.3959  -0.0582   6.2895  
# 
# Coefficients:
#   Estimate Std. Error   z value Pr(>|z|)    
# (Intercept)                 -4.194e+00  6.020e-02   -69.669  < 2e-16 ***
#   pm25                         5.460e-03  1.930e-04    28.283  < 2e-16 ***
#   nox                         -1.578e-04  1.453e-04    -1.086    0.277    
# poly(followupAge, 2)1        6.684e+02  2.288e+00   292.186  < 2e-16 ***
#   poly(followupAge, 2)2       -6.005e+02  2.415e+00  -248.680  < 2e-16 ***
#   as.factor(region)NE         -2.000e-02  1.339e-03   -14.937  < 2e-16 ***
#   as.factor(region)S           1.470e-01  1.522e-03    96.611  < 2e-16 ***
#   as.factor(region)W          -2.748e-01  2.324e-03  -118.281  < 2e-16 ***
#   as.factor(followupYear)2001 -7.635e-01  1.882e-03  -405.596  < 2e-16 ***
#   as.factor(followupYear)2002 -1.172e+00  1.872e-03  -626.042  < 2e-16 ***
#   as.factor(followupYear)2003 -1.458e+00  1.925e-03  -757.770  < 2e-16 ***
#   as.factor(followupYear)2004 -1.665e+00  1.960e-03  -849.463  < 2e-16 ***
#   as.factor(followupYear)2005 -1.893e+00  1.989e-03  -951.611  < 2e-16 ***
#   as.factor(followupYear)2006 -2.068e+00  2.091e-03  -988.738  < 2e-16 ***
#   as.factor(followupYear)2007 -2.205e+00  2.176e-03 -1013.243  < 2e-16 ***
#   as.factor(followupYear)2008 -2.339e+00  2.208e-03 -1059.318  < 2e-16 ***
#   as.factor(followupYear)2009 -2.442e+00  2.304e-03 -1060.166  < 2e-16 ***
#   as.factor(followupYear)2010 -2.539e+00  2.347e-03 -1081.582  < 2e-16 ***
#   as.factor(followupYear)2011 -2.692e+00  2.437e-03 -1104.606  < 2e-16 ***
#   as.factor(followupYear)2012 -3.296e+00  2.813e-03 -1171.726  < 2e-16 ***
#   as.factor(followupYear)2013 -3.478e+00  2.903e-03 -1198.017  < 2e-16 ***
#   as.factor(followupYear)2014 -3.629e+00  2.959e-03 -1226.490  < 2e-16 ***
#   as.factor(followupYear)2015 -3.304e+00  2.641e-03 -1251.107  < 2e-16 ***
#   as.factor(followupYear)2016 -5.506e+00  2.412e-03 -2282.828  < 2e-16 ***
#   as.factor(sex)2              2.627e-01  8.046e-04   326.509  < 2e-16 ***
#   as.factor(race)his          -8.562e-02  3.285e-03   -26.068  < 2e-16 ***
#   as.factor(race)oth          -3.229e-01  3.137e-03  -102.940  < 2e-16 ***
#   as.factor(race)wht          -1.254e-01  1.518e-03   -82.587  < 2e-16 ***
#   as.factor(dual)1             2.196e-01  1.077e-03   203.864  < 2e-16 ***
#   hispanic                    -7.524e-02  3.408e-03   -22.075  < 2e-16 ***
#   pct_blk                      5.245e-02  2.791e-03    18.791  < 2e-16 ***
#   bmi                          2.446e-02  5.015e-04    48.773  < 2e-16 ***
#   smoke                        1.490e-01  6.407e-03    23.250  < 2e-16 ***
#   income                       2.754e-06  4.142e-08    66.487  < 2e-16 ***
#   value                        2.256e-07  5.644e-09    39.978  < 2e-16 ***
#   poverty                      1.330e-02  9.172e-03     1.450    0.147    
# education                    5.578e-02  4.633e-03    12.041  < 2e-16 ***
#   popdensity                  -1.067e-06  5.740e-08   -18.593  < 2e-16 ***
#   ownerocc                    -2.068e-01  3.885e-03   -53.239  < 2e-16 ***
#   summer_tmmx                  1.019e-02  2.426e-04    42.012  < 2e-16 ***
#   summer_rmax                  6.015e-04  7.844e-05     7.667 1.75e-14 ***
#   winter_tmmx                 -3.999e-03  1.265e-04   -31.623  < 2e-16 ***
#   winter_rmax                 -4.606e-06  7.439e-05    -0.062    0.951    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 29996685  on 21143775  degrees of freedom
# Residual deviance: 13292303  on 21143733  degrees of freedom
# AIC: 25103703
# 
# Number of Fisher Scoring iterations: 6


########################## 2a. 5 yr age bins ###################################

ADRDdat5 <- ADRDdat[, nADRDEvents = sum(nADRDEvents)]


########################## 3. Smooth effects ###################################
############################ 3a. full cohort ###################################
m2 <- bam(nADRDEvents ~ 
            offset(log(I(nYears * N))) + # offset by total person-time at-risk
            s(pm25, bs = "cr", k = 3) +
            s(nox, bs = "cr", k = 3) +
            poly(followupAge, 2) +
            as.factor(region) + as.factor(followupYear) +
            as.factor(sex) + as.factor(race) + as.factor(dual) +
            hispanic + pct_blk +
            bmi + smoke + income + value + poverty +
            education + popdensity + ownerocc +
            summer_tmmx + summer_rmax + winter_tmmx + winter_rmax,
          data = ADRDdat,
          family = poisson(link = "log"),
          discrete = TRUE,
          nthreads = 16, 
          control = gam.control(trace = TRUE)) # print diagnostic output
summary(m2)
# Family: poisson 
# Link function: log 
# 
# Formula:
#   nADRDEvents ~ offset(log(I(nYears * N))) + s(pm25, bs = "cr", 
#                                                k = 3) + s(nox, bs = "cr", k = 3) + poly(followupAge, 2) + 
#   as.factor(region) + as.factor(followupYear) + as.factor(sex) + 
#   as.factor(race) + as.factor(dual) + hispanic + pct_blk + 
#   bmi + smoke + income + value + poverty + education + popdensity + 
#   ownerocc + summer_tmmx + summer_rmax + winter_tmmx + winter_rmax
# 
# Parametric coefficients:
#   Estimate Std. Error   z value Pr(>|z|)    
# (Intercept)                 -4.391e+00  6.128e-02   -71.649  < 2e-16 ***
#   poly(followupAge, 2)1        6.685e+02  2.287e+00   292.258  < 2e-16 ***
#   poly(followupAge, 2)2       -5.998e+02  2.415e+00  -248.363  < 2e-16 ***
#   as.factor(region)NE         -2.190e-02  1.339e-03   -16.348  < 2e-16 ***
#   as.factor(region)S           1.441e-01  1.523e-03    94.633  < 2e-16 ***
#   as.factor(region)W          -2.519e-01  2.395e-03  -105.162  < 2e-16 ***
#   as.factor(followupYear)2001 -7.618e-01  1.883e-03  -404.630  < 2e-16 ***
#   as.factor(followupYear)2002 -1.170e+00  1.874e-03  -624.334  < 2e-16 ***
#   as.factor(followupYear)2003 -1.457e+00  1.928e-03  -755.620  < 2e-16 ***
#   as.factor(followupYear)2004 -1.662e+00  1.967e-03  -844.936  < 2e-16 ***
#   as.factor(followupYear)2005 -1.893e+00  1.990e-03  -951.323  < 2e-16 ***
#   as.factor(followupYear)2006 -2.068e+00  2.097e-03  -986.211  < 2e-16 ***
#   as.factor(followupYear)2007 -2.205e+00  2.181e-03 -1010.935  < 2e-16 ***
#   as.factor(followupYear)2008 -2.337e+00  2.213e-03 -1055.655  < 2e-16 ***
#   as.factor(followupYear)2009 -2.434e+00  2.311e-03 -1053.283  < 2e-16 ***
#   as.factor(followupYear)2010 -2.539e+00  2.347e-03 -1081.787  < 2e-16 ***
#   as.factor(followupYear)2011 -2.693e+00  2.436e-03 -1105.238  < 2e-16 ***
#   as.factor(followupYear)2012 -3.299e+00  2.814e-03 -1172.164  < 2e-16 ***
#   as.factor(followupYear)2013 -3.474e+00  2.905e-03 -1195.963  < 2e-16 ***
#   as.factor(followupYear)2014 -3.626e+00  2.960e-03 -1224.875  < 2e-16 ***
#   as.factor(followupYear)2015 -3.303e+00  2.645e-03 -1248.895  < 2e-16 ***
#   as.factor(followupYear)2016 -5.507e+00  2.430e-03 -2266.701  < 2e-16 ***
#   as.factor(sex)2              2.626e-01  8.046e-04   326.310  < 2e-16 ***
#   as.factor(race)his          -8.527e-02  3.284e-03   -25.963  < 2e-16 ***
#   as.factor(race)oth          -3.250e-01  3.137e-03  -103.585  < 2e-16 ***
#   as.factor(race)wht          -1.251e-01  1.518e-03   -82.418  < 2e-16 ***
#   as.factor(dual)1             2.201e-01  1.077e-03   204.246  < 2e-16 ***
#   hispanic                    -7.640e-02  3.413e-03   -22.384  < 2e-16 ***
#   pct_blk                      4.704e-02  2.793e-03    16.843  < 2e-16 ***
#   bmi                          2.463e-02  5.020e-04    49.056  < 2e-16 ***
#   smoke                        1.481e-01  6.427e-03    23.046  < 2e-16 ***
#   income                       2.699e-06  4.146e-08    65.086  < 2e-16 ***
#   value                        2.136e-07  5.661e-09    37.742  < 2e-16 ***
#   poverty                      3.670e-02  9.180e-03     3.998 6.40e-05 ***
#   education                    4.265e-02  4.648e-03     9.177  < 2e-16 ***
#   popdensity                  -1.003e-06  5.752e-08   -17.435  < 2e-16 ***
#   ownerocc                    -1.958e-01  3.897e-03   -50.250  < 2e-16 ***
#   summer_tmmx                  1.052e-02  2.447e-04    42.992  < 2e-16 ***
#   summer_rmax                  7.195e-04  7.863e-05     9.149  < 2e-16 ***
#   winter_tmmx                 -3.822e-03  1.267e-04   -30.166  < 2e-16 ***
#   winter_rmax                 -3.374e-04  7.474e-05    -4.514 6.35e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(pm25) 1.831  1.972  596.9  <2e-16 ***
#   s(nox)  1.999  2.000 2323.0  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.406   Deviance explained = 35.1%
# fREML = 2.6075e+07  Scale est. = 1         n = 21143776
plot(m2, scale = 0, trans = exp)


####################### 3b. low-exposure ################################
m3 <- bam(nADRDEvents ~ 
            offset(log(I(nYears * N))) + # offset by total person-time at-risk
            s(pm25, bs = "cr", k = 3) +
            s(nox, bs = "cr", k = 3) +
            poly(followupAge, 2) +
            as.factor(region) + as.factor(followupYear) +
            as.factor(sex) + as.factor(race) + as.factor(dual) +
            hispanic + pct_blk +
            bmi + smoke + income + value + poverty +
            education + popdensity + ownerocc +
            summer_tmmx + summer_rmax + winter_tmmx + winter_rmax,
          data = ADRDdat, subset = which(ADRDdat$pm25 <= 12),
          family = poisson(link = "log"),
          discrete = TRUE,
          nthreads = 16, 
          control = gam.control(trace = TRUE)) # print diagnostic output
summary(m3)
# Family: poisson 
# Link function: log 
# 
# Formula:
#   nADRDEvents ~ offset(log(I(nYears * N))) + s(pm25, bs = "cr", 
#                                                k = 3) + s(nox, bs = "cr", k = 3) + s(followupAge, bs = "cr", 
#                                                                                      k = 3) + as.factor(region) + as.factor(followupYear) + as.factor(sex) + 
#   as.factor(race) + as.factor(dual) + hispanic + pct_blk + 
#   bmi + smoke + income + value + poverty + education + popdensity + 
#   ownerocc
# 
# Parametric coefficients:
#   Estimate Std. Error   z value Pr(>|z|)    
# (Intercept)                 -2.754e+00  1.764e-02  -156.131  < 2e-16 ***
#   as.factor(region)NE          2.133e-02  1.637e-03    13.031  < 2e-16 ***
#   as.factor(region)S           1.870e-01  1.391e-03   134.408  < 2e-16 ***
#   as.factor(region)W          -2.435e-01  2.147e-03  -113.408  < 2e-16 ***
#   as.factor(followupYear)2001 -7.470e-01  3.076e-03  -242.828  < 2e-16 ***
#   as.factor(followupYear)2002 -1.160e+00  3.023e-03  -383.835  < 2e-16 ***
#   as.factor(followupYear)2003 -1.443e+00  3.096e-03  -466.016  < 2e-16 ***
#   as.factor(followupYear)2004 -1.648e+00  3.012e-03  -547.309  < 2e-16 ***
#   as.factor(followupYear)2005 -1.884e+00  3.242e-03  -581.084  < 2e-16 ***
#   as.factor(followupYear)2006 -2.046e+00  3.050e-03  -670.717  < 2e-16 ***
#   as.factor(followupYear)2007 -2.189e+00  3.209e-03  -682.070  < 2e-16 ***
#   as.factor(followupYear)2008 -2.297e+00  2.923e-03  -785.785  < 2e-16 ***
#   as.factor(followupYear)2009 -2.407e+00  2.956e-03  -814.248  < 2e-16 ***
#   as.factor(followupYear)2010 -2.490e+00  2.959e-03  -841.626  < 2e-16 ***
#   as.factor(followupYear)2011 -2.628e+00  2.944e-03  -892.872  < 2e-16 ***
#   as.factor(followupYear)2012 -3.251e+00  3.347e-03  -971.337  < 2e-16 ***
#   as.factor(followupYear)2013 -3.434e+00  3.434e-03  -999.871  < 2e-16 ***
#   as.factor(followupYear)2014 -3.595e+00  3.483e-03 -1032.185  < 2e-16 ***
#   as.factor(followupYear)2015 -3.262e+00  3.214e-03 -1014.976  < 2e-16 ***
#   as.factor(followupYear)2016 -5.476e+00  3.070e-03 -1783.361  < 2e-16 ***
#   as.factor(sex)2              2.544e-01  9.967e-04   255.298  < 2e-16 ***
#   as.factor(race)his          -7.100e-02  3.959e-03   -17.935  < 2e-16 ***
#   as.factor(race)oth          -3.091e-01  3.860e-03   -80.076  < 2e-16 ***
#   as.factor(race)wht          -1.097e-01  2.042e-03   -53.700  < 2e-16 ***
#   as.factor(dual)1             2.482e-01  1.415e-03   175.362  < 2e-16 ***
#   hispanic                    -6.685e-02  3.811e-03   -17.541  < 2e-16 ***
#   pct_blk                      7.686e-02  3.835e-03    20.039  < 2e-16 ***
#   bmi                          3.397e-02  6.133e-04    55.391  < 2e-16 ***
#   smoke                        1.131e-01  7.839e-03    14.428  < 2e-16 ***
#   income                       3.103e-06  4.790e-08    64.781  < 2e-16 ***
#   value                        1.233e-07  6.193e-09    19.907  < 2e-16 ***
#   poverty                      1.632e-02  1.151e-02     1.418    0.156    
# education                    4.124e-02  5.928e-03     6.957 3.47e-12 ***
#   popdensity                  -1.501e-06  8.312e-08   -18.061  < 2e-16 ***
#   ownerocc                    -2.417e-01  4.782e-03   -50.549  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df  Chi.sq p-value    
# s(pm25)        1.998      2   688.8  <2e-16 ***
#   s(nox)         1.999      2  1205.4  <2e-16 ***
#   s(followupAge) 2.000      2 66680.3  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.367   Deviance explained = 33.8%
# fREML = 1.9709e+07  Scale est. = 1         n = 16239121
plot(m3, scale = 0, trans = exp)




####################### 3c. sex stratification ################################
m4 <- bam(nADRDEvents ~ 
            offset(log(I(nYears * N))) + # offset by total person-time at-risk
            s(pm25, bs = "cr", k = 3, by = as.factor(sex)) +
            s(nox, bs = "cr", k = 3) +
            poly(followupAge, 2) +
            as.factor(region) + as.factor(followupYear) +
            as.factor(sex) + as.factor(race) + as.factor(dual) +
            hispanic + pct_blk +
            bmi + smoke + income + value + poverty +
            education + popdensity + ownerocc +
            summer_tmmx + summer_rmax + winter_tmmx + winter_rmax,
          data = ADRDdat,
          family = poisson(link = "log"),
          discrete = TRUE,
          nthreads = 16, 
          control = gam.control(trace = TRUE)) # print diagnostic output
summary(m4)
# Family: poisson 
# Link function: log 
# 
# Formula:
#   nADRDEvents ~ offset(log(I(nYears * N))) + s(pm25, bs = "cr", 
#                                                k = 3, by = as.factor(sex)) + s(nox, bs = "cr", k = 3) + 
#   poly(followupAge, 2) + as.factor(region) + as.factor(followupYear) + 
#   as.factor(sex) + as.factor(race) + as.factor(dual) + hispanic + 
#   pct_blk + bmi + smoke + income + value + poverty + education + 
#   popdensity + ownerocc + summer_tmmx + summer_rmax + winter_tmmx + 
#   winter_rmax
# 
# Parametric coefficients:
#   Estimate Std. Error   z value Pr(>|z|)    
# (Intercept)                 -4.406e+00  6.104e-02   -72.184  < 2e-16 ***
#   poly(followupAge, 2)1        6.674e+02  2.288e+00   291.707  < 2e-16 ***
#   poly(followupAge, 2)2       -5.989e+02  2.415e+00  -247.969  < 2e-16 ***
#   as.factor(region)NE         -2.177e-02  1.339e-03   -16.257  < 2e-16 ***
#   as.factor(region)S           1.440e-01  1.522e-03    94.588  < 2e-16 ***
#   as.factor(region)W          -2.519e-01  2.386e-03  -105.608  < 2e-16 ***
#   as.factor(followupYear)2001 -7.617e-01  1.883e-03  -404.598  < 2e-16 ***
#   as.factor(followupYear)2002 -1.170e+00  1.873e-03  -624.447  < 2e-16 ***
#   as.factor(followupYear)2003 -1.456e+00  1.926e-03  -756.000  < 2e-16 ***
#   as.factor(followupYear)2004 -1.662e+00  1.965e-03  -845.976  < 2e-16 ***
#   as.factor(followupYear)2005 -1.893e+00  1.990e-03  -951.536  < 2e-16 ***
#   as.factor(followupYear)2006 -2.067e+00  2.094e-03  -987.109  < 2e-16 ***
#   as.factor(followupYear)2007 -2.204e+00  2.179e-03 -1011.640  < 2e-16 ***
#   as.factor(followupYear)2008 -2.336e+00  2.212e-03 -1056.331  < 2e-16 ***
#   as.factor(followupYear)2009 -2.434e+00  2.311e-03 -1053.403  < 2e-16 ***
#   as.factor(followupYear)2010 -2.539e+00  2.347e-03 -1081.639  < 2e-16 ***
#   as.factor(followupYear)2011 -2.693e+00  2.436e-03 -1105.113  < 2e-16 ***
#   as.factor(followupYear)2012 -3.298e+00  2.814e-03 -1172.245  < 2e-16 ***
#   as.factor(followupYear)2013 -3.474e+00  2.905e-03 -1196.119  < 2e-16 ***
#   as.factor(followupYear)2014 -3.626e+00  2.960e-03 -1224.884  < 2e-16 ***
#   as.factor(followupYear)2015 -3.303e+00  2.643e-03 -1249.564  < 2e-16 ***
#   as.factor(followupYear)2016 -5.507e+00  2.422e-03 -2273.780  < 2e-16 ***
#   as.factor(sex)2              2.683e-01  2.878e-03    93.217  < 2e-16 ***
#   as.factor(race)his          -8.492e-02  3.284e-03   -25.856  < 2e-16 ***
#   as.factor(race)oth          -3.245e-01  3.137e-03  -103.439  < 2e-16 ***
#   as.factor(race)wht          -1.252e-01  1.518e-03   -82.462  < 2e-16 ***
#   as.factor(dual)1             2.193e-01  1.078e-03   203.469  < 2e-16 ***
#   hispanic                    -7.620e-02  3.412e-03   -22.336  < 2e-16 ***
#   pct_blk                      4.691e-02  2.793e-03    16.796  < 2e-16 ***
#   bmi                          2.467e-02  5.019e-04    49.149  < 2e-16 ***
#   smoke                        1.473e-01  6.420e-03    22.947  < 2e-16 ***
#   income                       2.701e-06  4.146e-08    65.158  < 2e-16 ***
#   value                        2.139e-07  5.657e-09    37.803  < 2e-16 ***
#   poverty                      3.628e-02  9.179e-03     3.953 7.73e-05 ***
#   education                    4.336e-02  4.645e-03     9.335  < 2e-16 ***
#   popdensity                  -1.008e-06  5.751e-08   -17.525  < 2e-16 ***
#   ownerocc                    -1.964e-01  3.895e-03   -50.426  < 2e-16 ***
#   summer_tmmx                  1.056e-02  2.438e-04    43.303  < 2e-16 ***
#   summer_rmax                  7.342e-04  7.859e-05     9.342  < 2e-16 ***
#   winter_tmmx                 -3.836e-03  1.266e-04   -30.298  < 2e-16 ***
#   winter_rmax                 -3.369e-04  7.474e-05    -4.508 6.54e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df  Chi.sq  p-value    
#   s(pm25):as.factor(sex)1 1.952  1.997   53.37 1.52e-12 ***
#   s(pm25):as.factor(sex)2 1.292  1.498  357.49  < 2e-16 ***
#   s(nox)                  2.000  2.000 2345.03  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.406   Deviance explained = 35.1%
# fREML = 2.6075e+07  Scale est. = 1         n = 21143776
plot(m4, scale = 0, trans = exp)



# PM2.5 plot
pm25_plot <- plot(m2, select = 1, xlim = quantile(ADRDdat$pm25, c(0.005, 0.995), na.rm = T))
plot(m2, select = 1, trans = exp, scale = 0, shift = -(pm25_plot[[1]]$fit[1]),
     xlim = quantile(ADRDdat$pm25, c(0.005, 0.995), na.rm = T))
abline(v = quantile(ADRDdat$pm25, c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = T))
abline(h = exp(pm25_plot[[1]]$fit[sapply(c(0.1, 0.25, 0.5, 0.75, 0.9), function(q) {
  which.min(abs(pm25_plot[[1]]$x - quantile(ADRDdat$pm25, q, na.rm = T)))
})] - (pm25_plot[[1]]$fit[1])))

# NOx plot
nox_plot <- plot(m2, select = 2, xlim = quantile(ADRDdat$nox, c(0.005, 0.995), na.rm = T))
plot(m2, select = 2, trans = exp, scale = 0, shift = -(nox_plot[[2]]$fit[1]),
     xlim = quantile(ADRDdat$nox, c(0.005, 0.995), na.rm = T))
abline(v = quantile(ADRDdat$nox, c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = T))
abline(h = exp(nox_plot[[2]]$fit[sapply(c(0.1, 0.25, 0.5, 0.75, 0.9), function(q) {
  which.min(abs(nox_plot[[2]]$x - quantile(ADRDdat$nox, q, na.rm = T)))
})] - (nox_plot[[2]]$fit[1])))

# Smooth age plot
plot(m2, select = 3, scale = 0, trans = exp)
