###############################################################################
# Project: Causal effect of air pollution on first ADRD hospitalization       #
# Code: get ADRD study population, exclude problematic ID                     #
# Input: "no_crosswalk_no_death_ids.fst"                                      #
# Input: "ADRD`type`_`combined`.fst"                                          #
# Input: "no_crosswalk_no_death_ids.fst"                                      #
# Output: "First_hosp_yr.fst" - IDs and first year of ADRD admission          #
# Author: Daniel Mork                                                         #
# Date: 2021-12-21                                                            #
###############################################################################
rm(list = ls())
gc()
############################# 0. Setup ########################################
library(data.table)
library(fst)
library(NSAPHutils)

setDTthreads(threads = 4)

dir_data <- "/nfs/home/D/dam9096/shared_space/ci3_analysis/dmork/Data/Causal_ADRD/"
dir_input_crosswalk <- "/nfs/home/D/dam9096/shared_space/ci3_health_data/medicare/id_crosswalk/"

################ 1. combine hospitalization files together #####################
ADRDhosp <- rbindlist(list(read_fst(paste0(dir_data, "ADRDprimary_combined.fst"), as.data.table = TRUE),
                           read_fst(paste0(dir_data, "ADRDsecondary_combined.fst"), as.data.table = TRUE))
names(ADRDhosp)
# [1] "QID"            "ADATE"          "DDATE"          "DIAG1"         
# [5] "DIAG2"          "DIAG3"          "DIAG4"          "DIAG5"         
# [9] "DIAG6"          "DIAG7"          "DIAG8"          "DIAG9"         
# [13] "DIAG10"         "year"           "ADRD_primary"   "ADRD_secondary"
dim(ADRDhosp) 
# [1] 19502539       16

dup <- duplicated(ADRDhosp)
sum(dup) 
# [1] 1540072
ADRDhosp <- unique(ADRDhosp)
dim(ADRDhosp) 
# [1] 17962467       16
ADRDhosp[, .N, by = c("year", "ADRD_primary")]
# year ADRD_primary       N
# 1: 2000         TRUE   95268
# 2: 2001         TRUE  101102
# 3: 2002         TRUE  102562
# 4: 2003         TRUE  105139
# 5: 2004         TRUE  107191
# 6: 2005         TRUE  108995
# 7: 2006         TRUE  101345
# 8: 2007         TRUE   98194
# 9: 2008         TRUE   93530
# 10: 2009         TRUE   89044
# 11: 2010         TRUE   89840
# 12: 2011         TRUE   83264
# 13: 2012         TRUE   81043
# 14: 2013         TRUE   75812
# 15: 2014         TRUE   73391
# 16: 2015         TRUE   72511
# 17: 2016         TRUE   61841
# 18: 2000        FALSE  863221
# 19: 2001        FALSE  922879
# 20: 2002        FALSE  973820
# 21: 2003        FALSE 1005547
# 22: 2004        FALSE 1001763
# 23: 2005        FALSE 1008598
# 24: 2006        FALSE  979714
# 25: 2007        FALSE  981911
# 26: 2008        FALSE  996064
# 27: 2009        FALSE  956465
# 28: 2010        FALSE  989965
# 29: 2011        FALSE 1034065
# 30: 2012        FALSE  974102
# 31: 2013        FALSE  935912
# 32: 2014        FALSE  915231
# 33: 2015        FALSE  964446
# 34: 2016        FALSE  918692
# year ADRD_primary       N

ADRDhosp <- ADRDhosp[, .(QID, ADATE, year)]
ADRDhosp <- unique(ADRDhosp)
ADRDhosp[, ADATE:=NULL]
ADRDhosp <- unique(ADRDhosp)
dim(ADRDhosp)
# [1] 13132546        2

ADRDhosp_first <- ADRDhosp[, list(firstADRDyr = min(year)), by = .(QID)]
setDT(ADRDhosp_first)
dim(ADRDhosp_first)
# [1] 8325790       2

ADRDhosp_first[, .N, by = "firstADRDyr"]
# firstADRDyr      N
# 1:        2000 707931
# 2:        2001 594763
# 3:        2002 556658
# 4:        2003 534299
# 5:        2004 511654
# 6:        2005 498213
# 7:        2006 474277
# 8:        2007 469835
# 9:        2008 459689
# 10:        2009 438021
# 11:        2010 455192
# 12:        2011 465374
# 13:        2012 439492
# 14:        2013 428435
# 15:        2014 424992
# 16:        2015 449241
# 17:        2016 417724

######################### 2. exclude problematic IDs ##########################
probIDs <- read_fst(paste0(dir_input_crosswalk, "no_crosswalk_no_death_ids.fst"), 
                    as.data.table = T)
ADRDhosp_first <- ADRDhosp_first[!(QID %in% probIDs[, old_id]),] # exclude problematic IDs
dim(ADRDhosp_first)
# [1] 8323618       2

############## 3. save first hospitalization INFO ##############################
ADRDhosp_first$QID <- factor(ADRDhosp_first$QID)
write_fst(ADRDhosp_first, paste0(dir_data, "First_hosp_yr.fst"))

