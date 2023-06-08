##### Setup #####

library(data.table)
library(fst)
library(mgcv)

# get directories, classifications of variables, and formulas
dir_code <- "~/nsaph_projects/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/git/code/"
source(paste0(dir_code, "constants.R"))

# set exposure
exposure_name <- "no2"

# set parameters for this computing job
n_cores <- 16

# get data
zip_year_data_with_strata <- read_fst(paste0(dir_data, "analysis/",
                                             exposure_name, "/",
                                             "zip_year_data_with_strata_trimmed_0.05_0.95.fst"),
                                      as.data.table = T)

# make sure categorical variables are factors
zip_year_data_with_strata[, `:=`(zip = as.factor(zip),
                                 year = as.factor(year),
                                 age_grp = as.factor(age_grp),
                                 sex = as.factor(sex),
                                 race = as.factor(race),
                                 dual = as.factor(dual))]

# run model (Poisson regression)
cl <- parallel::makeCluster(n_cores, type = "PSOCK")
bam_associational <- bam(as.formula(paste("Y ~", paste(c("w",
                                                         strata_vars,
                                                         zip_var_names),
                                                       collapse = "+", sep = ""))),
                         data = zip_year_data_with_strata,
                         offset = log(n_persons * n_years),
                         family = poisson(link = "log"),
                         samfrac = 0.05,
                         chunk.size = 5000,
                         control = gam.control(trace = TRUE),
                         nthreads = n_cores,
                         cluster = cl)
parallel::stopCluster(cl)
cat(paste(exposure_name, "Poisson Regression", bam_associational$coefficients["w"], sep = ","),
    sep = "\n",
    file = paste0(dir_results, "parametric_results/coef_for_exposure.txt"),
    append = TRUE)


### Sensitivity analysis: thin-plate spline outcome model ###

# fit model
cl <- parallel::makeCluster(n_cores, type = "PSOCK")
bam_smooth_associational <- 
  bam(as.formula(paste("Y ~", paste(c("s(w, bs = 'cr', k = 4, fx = TRUE)",
                                      strata_vars,
                                      zip_var_names),
                                    collapse = "+", sep = ""))),
                         data = zip_year_data_with_strata,
                         offset = log(n_persons * n_years),
                         family = poisson(link = "log"),
                         samfrac = 0.05,
                         chunk.size = 5000,
                         control = gam.control(trace = TRUE),
                         nthreads = n_cores,
                         cluster = cl)
parallel::stopCluster(cl)

# estimate counterfactual for every year-zip-strata, calculate ATE
potential_data <- copy(zip_year_data_with_strata)
data_prediction <- 
  rbindlist(lapply(seq(min(zip_year_data_with_strata$w), 
                       max(zip_year_data_with_strata$w), 
                       length.out = 100), function(pot_exp) {
                         
     # Get potential data if all had same potential exposure
     potential_data[, w := pot_exp]
     return(data.table(name = exposure_name,
                       w = pot_exp,
                       ate = mean(predict(bam_smooth_associational, 
                                          newdata = potential_data, 
                                          type = "response"))))
   }))
plot(I(1e5*ate)~w,data_prediction, type = 'l')
save(data_prediction, file = paste0(dir_results, exposure_name, "_assoc_smooth.rda"))