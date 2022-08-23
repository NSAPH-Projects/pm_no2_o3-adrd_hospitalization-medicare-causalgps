###############################################################################
# Project: Causal effect of air pollution on first ADRD hospitalization       #
# Code: visualize AD/ADRD hospitalization data by histogram                   #
# Input: hospitalization files                                                #
# Author: Daniel Mork                                                         #
# Date: 2022-01-27                                                            #
###############################################################################
rm(list = ls())
gc()
cores = 2
############################# 0. Setup ########################################
library(data.table)
library(fst)
library(NSAPHutils)
library(lubridate)
library(icd)

setDTthreads(threads = cores)
dir_hospital <- "/nfs/home/D/dam9096/shared_space/ci3_health_data/medicare/gen_admission/1999_2016/targeted_conditions/cache_data/admissions_by_year/"
dir_data <- "/nfs/home/D/dam9096/shared_space/ci3_analysis/dmork/Data/Causal_ADRD/"

# Histogram of hospitalizations by year
anyAD <- read_fst(paste0(dir_data, "Hosp_outcomes/AD_any.fst"), as.data.table = TRUE)
anyADRD <- read_fst(paste0(dir_data, "Hosp_outcomes/ADRD_any.fst"), as.data.table = TRUE)

# Primary or secondary AD/ADRD hosp
hist(anyAD[, .(year = first(year)), by = c("QID", "ADATE")][, year], xlab = "Year", main = "AD Hospitalizations")
hist(anyADRD[, .(year = first(year)), by = c("QID", "ADATE")][, year], xlab = "Year", main = "ADRD Hospitalizations")

# Primary only AD/ADRD hosp
hist(anyAD[DIAG=="DIAG1", year, by = c("QID", "ADATE")][, year], xlab = "Year", 
     main = "AD Hospitalizations")
hist(anyADRD[DIAG=="DIAG1", year, by = c("QID", "ADATE")][, year], xlab = "Year", 
     main = "ADRD Hospitalizations")

# By ICD code, primary or secondary
ICD <- unique(anyADRD$ICD)
for(c in ICD) {
  hist(anyADRD[ICD==c, .(year = first(year)), by = c("QID", "ADATE")][, year], xlab = "Year", 
       main = c)
}
