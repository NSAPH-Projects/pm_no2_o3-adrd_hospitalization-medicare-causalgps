###############################################################################
# Project: Causal effect of air pollution on first ADRD hospitalization       #
# Code: extract ADRD hopspitalization info from hospitalization files         #
# Input: hospitalization files                                                #
# Output: "ADRD`type`_`year`.fst"                                             #
# Output: "ADRD`type`_`combined`.fst"                                         #
# Author: Daniel Mork                                                         #
# Date: 2021-12-18                                                            #
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

######################## 1. ICD code info #####################################
# using lists from CCW: https://www2.ccwdata.org/web/guest/condition-categories
outcomes <- list()
outcomes[["AD"]] <- list(
  "icd9" = c(children("3310")),
  "icd10" = c(children("G300"), children("G301"), children("G308"), children("G309"))
)
outcomes[["ADRD"]] <- list(
  "icd9" = c(children("2900"), children("2901"), children("2902"), children("2903"), children("2904"),
             children("2941"), children("2942"), children("2948"), children("797"),
             children("3310"), children("3311"), children("3312"), children("3317")),
  "icd10" = c(children("F015"), children("F028"), children("F039"), children("F04"),
              children("G138"), children("F05"), children("F061"), children("F068"),
              outcomes[["AD"]][["icd10"]], children("G311"), children("G312"),
              children("G310"), children("G94"), children("R4184"), children("R54"))
)
outcomeDT <- rbindlist(list(
  data.table(icd = as.character(outcomes[["AD"]][["icd9"]]), AD = TRUE, RD = FALSE),
  data.table(icd = as.character(outcomes[["AD"]][["icd10"]]), AD = TRUE, RD = FALSE),
  data.table(icd = as.character(outcomes[["ADRD"]][["icd9"]]), AD = FALSE, RD = TRUE),
  data.table(icd = as.character(outcomes[["ADRD"]][["icd10"]]), AD = FALSE, RD = TRUE)))
outcomeDT <- unique(outcomeDT[, .(AD = any(AD), ADRD = AD | RD), by = icd])


# From: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8075520/
# mild cognitive impairment (ICD-9: 331.81, 294.9; ICD-10: G31.83, F09), 
# Alzheimer's disease (ICD-9: 331.0; ICD-10: G30.0, G30.1, G30.8, G30.9), 
# vascular dementia (ICD-9: 290.40, 290.41; ICD-10: F01.50, F01.51), 
# Lewy body dementia (ICD-9: 331.82; ICD-10: G31.83), 
# Frontotemporal Dementia (ICD-9: 331.19; ICD-10: G31.09), and 
# primary progressive aphasia (ICD-9: 331.11; ICD-10: G31.01)

# Final list
# outcomes[["ADRD"]][["icd9"]] <- c(children("290"), # Dementias
#                                   children("2941"), # Dementia, classified elsewhere
#                                   children("2942"), # Dementia, unspecified
#                                   children("2948"), # Other persistent mental disorder
#                                   children("2949"), # Unspecified persistent mental disorder
#                                   children("331"), # Alzheimer's disease
#                                   children("797"))
# outcomes[["ADRD"]][["icd10"]] <- c(children("F01"), # Vascular dementia
#                                    children("F02"), # Dementia in other diseases
#                                    children("F03"), # Unspecified dementia 
#                                    children("F09"), # MCI
#                                    children("G30"), # Alzheimer's disease
#                                    children("G310"), # Frontotemoral dementia
#                                    children("G311"), # Senile dementia
#                                    "G3183", # Dementia with Lewy bodies
#                                    "G3184") # MCI, so stated



##################### 2. extract hospitalization info #########################
## uncomment next line to clear out old data in case of re-run
# file.remove(list.files(dir_data, pattern = ".fst", full.names = T))
ADRDhosp <- data.table()
for (year_ in 2000:2016) {
  cat("Loading", year_, "hospitalization file... \n")
  admissions <- read_fst(paste0(dir_hospital, "admissions_", year_, ".fst"),
                         c("QID", "ADATE", "DDATE", paste0("DIAG", 1:10)),
                         as.data.table = TRUE)
  admissions[, ADATE := dmy(ADATE)][, DDATE := dmy(DDATE)]
  admissions[, year := year(ADATE)]
  admissions <- admissions[year == year_]
  admissions <- melt(admissions, measure.vars = paste0("DIAG", 1:10),
                     variable.name = "DIAG", value.name = "ICD")
  admissions[, primary := (DIAG == "DIAG1")]
  admissions <- merge(admissions, outcomeDT,
                      by.x = "ICD", by.y = "icd", all.x = TRUE)[AD | ADRD]
  print(admissions[, .N, by = c("primary", "AD", "ADRD")])
  ADRDhosp <- rbindlist(list(ADRDhosp, admissions))
  rm(admissions); gc()
} # end open yearly hosp. files

########################### 3. save hospitalization data #######################
# filter primary/secondary AD/ADRD hospitalizations
primaryAD <- ADRDhosp[primary & AD]
write_fst(primaryAD, paste0(dir_data, "Hosp_outcomes/AD_primary.fst"))
primaryADRD <- ADRDhosp[primary & ADRD]
write_fst(primaryADRD, paste0(dir_data, "Hosp_outcomes/ADRD_primary.fst"))
# these may contain 'duplicate' rows because of multiple AD/ADRD diag codes
anyAD <- ADRDhosp[AD == TRUE]
write_fst(anyAD, paste0(dir_data, "Hosp_outcomes/AD_any.fst"))
anyADRD <- ADRDhosp[ADRD == TRUE]
write_fst(anyADRD, paste0(dir_data, "Hosp_outcomes/ADRD_any.fst"))

# filter to first hospitalization only
firstPrimaryAD <- primaryAD[, .(ADATE = min(ADATE), year = first(year)), by = QID]
write_fst(firstPrimaryAD, paste0(dir_data, "Hosp_outcomes/First_hosp_AD_primary.fst"))
firstPrimaryADRD <- primaryADRD[, .(ADATE = min(ADATE), year = first(year)), by = QID]
write_fst(firstPrimaryADRD, paste0(dir_data, "Hosp_outcomes/First_hosp_ADRD_primary.fst"))
firstAnyAD <- anyAD[, .(ADATE = min(ADATE), year = first(year)), by = QID]
write_fst(firstAnyAD, paste0(dir_data, "Hosp_outcomes/First_hosp_AD_any.fst"))
firstAnyADRD <- anyADRD[, .(ADATE = min(ADATE), year = first(year)), by = QID]
write_fst(firstAnyADRD, paste0(dir_data, "Hosp_outcomes/First_hosp_ADRD_any.fst"))