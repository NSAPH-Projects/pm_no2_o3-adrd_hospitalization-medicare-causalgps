rm(list = ls())
gc()

##### 0. Setup #####
library(data.table)
library(fst)
# library(NSAPHutils)
library(CausalGPS)

setDTthreads(threads = 16)

dir_data <- "/nfs/home/M/miq336/shared_space/ci3_analysis/dmork/Data/Causal_ADRD/"
ADRD_agg <- read_fst(paste0(dir_data, "analysis/ADRD_complete.fst"), as.data.table = TRUE)

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

ADRD_agg <- ADRD_agg[ADRD_year - ffs_entry_year >= 2, ]

# generate_pseudo_pop
# estimate_gps
# estimate_erf
# generate_syn_data

dt_subset <- ADRD_agg[sexM == F & race_cat == "black" & any_dual == F & ADRD_age == 70]
Y <- dt_subset$n_ADRDhosp / dt_subset$n_persons / dt_subset$n_years
w <- dt_subset$pm25
c <- dt_subset
c <- c[, `:=`(n_ADRDhosp = NULL, n_persons = NULL, n_years = NULL, pm25 = NULL)]

data_with_gps <- estimate_gps(Y,
                              w,
                              c,
                              pred_model = "sl",
                              gps_model = "parametric",
                              internal_use = FALSE,
                              params = list(xgb_max_depth = c(3,4,5),
                                            xgb_rounds = c(10,20,30,40)),
                              nthread = 1,                                
                              sl_lib = c("m_xgboost")
)
