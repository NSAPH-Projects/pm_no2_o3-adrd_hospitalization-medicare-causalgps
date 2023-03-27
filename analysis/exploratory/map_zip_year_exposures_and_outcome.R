# Memory required to run: More than 16 GB?

# Loading packages
library(data.table)
library(fst)
library(sf)
library(dplyr)

# directories for data, code, and results
dir_data <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/data/"
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/code/"
dir_results <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/results/"

# directory for ZIP code shapefiles
dir_shape <- "/n/dominici_nsaph_l3/Lab/data/shapefiles/zip_shape_files/Zipcode_Info/polygon/"

# read in full data
ADRD_agg_lagged <- read_fst(paste0(dir_data, "analysis/ADRD_complete_tv.fst"),
                            as.data.table = TRUE)
ADRD_agg_lagged[, zip := formatC(zip, width = 5, flag = "0")] # turn ZIP code into character and make sure there are leading zeroes

# read in ZIP code shapefiles, as a list by year
zip_shape_list <- lapply(0:15, function(i) read_sf(paste0(dir_shape, "ESRI", formatC(i, width = 2, flag = "0"), "USZIP5_POLY_WGS84.shp")))

# get exposure data as 1 observation per ZIP year; note that exposure year is one year behind patient year (long-term exposure)
zip_year_expos <- subset(ADRD_agg_lagged, select = c("pm25", "no2", "ozone_summer",
                                                     "zip", "year"))
zip_year_expos <- unique(zip_year_expos, by = c("zip", "year"))
zip_year_expos[, year := year - 1]

# split exposure data into list by year (2000-2015)
zip_expos_list <- split(zip_year_expos, by = "year")

# get number of patients (= number of person years since it's 1 year) and events in 2016 by ZIP code
# note: sum the patients and events to unstratify them
zip_patients_2016 <- subset(ADRD_agg_lagged, subset = ADRD_agg_lagged$year == "2016", select = c("n_hosp", "n_persons", "zip"))
zip_patients_2016 <- zip_patients_2016[, .(n_hosp = sum(n_hosp), n_persons = sum(n_persons)), by = .(zip)]
zip_patients_2016[, hosp_rate := n_hosp / n_persons]

# get patients' total outcome data: sum number of person-years and number of ADRD hospitalizations in each ZIP code (over all years)
zip_year_patients <- subset(ADRD_agg_lagged, select = c("n_hosp", "n_persons",
                                                        "zip", "year"))
zip_year_patients <- zip_year_patients[, .(n_hosp = sum(n_hosp), n_person_years = sum(n_persons)), by = .(zip)]
zip_year_patients[, overall_rate := n_hosp / n_person_years] # total rate in each ZIP code over 2001-2016

# merge exposure data with ZIP code shapefiles by year
zip_expos_merged <- lapply(0:15, function(i) zip_shape_list[[i + 1]] %>%
                             dplyr::right_join(zip_expos_list[[i + 1]], by = c("ZIP" = "zip")))

# merge patient data with 2016 ZIP code shapefile
zip_patients_2016_merged <- zip_shape_list[[16]] %>%
  dplyr::right_join(zip_patients_2016, by = c("ZIP" = "zip"))
zip_patients_merged <- zip_shape_list[[16]] %>%
  dplyr::right_join(zip_year_patients, by = c("ZIP" = "zip"))

# plot exposures in 2000
plot(zip_expos_merged[[1]]["pm25"], main = "PM2.5 (micrograms/m^3) in 2000")
plot(zip_expos_merged[[1]]["no2"], main = "NO2 (ppb) in 2000")
plot(zip_expos_merged[[1]]["ozone_summer"], main = "Summer ozone (ppb) in 2000")

# plot exposures in 2015
plot(zip_expos_merged[[16]]["pm25"], main = "PM2.5 (micrograms/m^3) in 2015")
plot(zip_expos_merged[[16]]["no2"], main = "NO2 (ppb) in 2015")
plot(zip_expos_merged[[16]]["ozone_summer"], main = "Summer ozone (ppb) in 2015")

# plot number of patients and number of ADRD hospitalizations in 2016
plot(zip_patients_2016_merged["n_persons"], main = "Number of patients in 2016")
plot(zip_patients_2016_merged["n_hosp"], main = "Number of hospitalizations with ADRD in 2016")
plot(zip_patients_2016_merged["hosp_rate"], main = "Rate of hospitalizations with ADRD in 2016")

# plot total (sum of) person-years and ADRD hospitalizations in 2001-2016
plot(zip_patients_merged["n_person_years"], main = "Total number of person-years over 2001-2016")
plot(zip_patients_merged["n_hosp"], main = "Total number of hospitalizations with ADRD over 2001-2016")
plot(zip_patients_merged["overall_rate"], main = "Overall rate of hospitalizations with ADRD over 2001-2016")