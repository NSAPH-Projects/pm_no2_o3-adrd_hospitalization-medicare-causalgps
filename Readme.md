**Project**: Causal_ADRD

**Authors**: Daniel Mork, Michelle Qin, Naeem Khoshnevis

**Created**: Jan. 21, 2022

**Overview**: Causal inference analysis of long-term air pollution (PM2.5, NO2, summer ozone) effect on first hospitalization with Alzheimer's disease and related dementias (ADRD).

**Data**: All Medicare beneficiaries continuously enrolled in FFS plan between 2001 and 2016, until first ADRD hospitalization or change to HMO plan or death, aggregated by entry year, event year, ZIP code, 5-year-binned age group, race, sex, and Medicaid eligibility.

**Description of Files and Folders**: 

`code/`

  - *constants.R*: Hard-codes directory paths, covariate names, and formulas for outcome models.

`code/aggregation/`   

These six files create the three datasets, one for each exposure, to be analyzed.   

  - *1. Medicare FFS enrollment.R*: Reads denominator files, gets 2000-2016 FFS data, defines variables like race and sex as at entry into FFS.
  - *2. ADRD hospitalization data.R*: Gets first ADRD hospitalization year, age, and ZIP code; merge into qid entry exit data.
  - *3. Aggregate data sets.R*: Aggregates individuals by FFS entry year (called `cohort`), `zip`, `age_grp`, `year`, `sex`, `race`, and Medicaid eligibility (called `dual`). Each observation is now an "aggregated patient unit".
  - *4. Year-zip confounders.R:* Gets confounders and exposures by zip and year.
  - *5. Merge and clean.R*: Merges patient data with confounders from current year and exposures from previous year, using complete cases.
  - *6. Set up trimmed and ZIP data.R*: For each of PM2.5, NO2, summer ozone at the ZIP code and annual level, trims the dataset at the exposure's 5th and 95th percentiles, keeping the middle 90\% of aggregated patient units.
    
`code/analysis/`
  
  - *helper_functions.R*: Defines some useful functions for transforming variables, checking covariate balance, causal inference methods, etc.
  
  - *setup_parametric_results_tables.R*: Create a .txt file to collect regression results and a .txt file to collect m-out-of-n bootstrap standard errors.
  
  - *associational_model.R*: Poisson regression model, used as a baseline to compare causal inference results.
  
  - *gps_weighting_by_zip_year.R*: GPS weighting analysis.
  
  - *gps_matching_by_zip_year.R*: GPS matching analysis.
  
`code/analysis/bootstrap_code/`
  
  - *setup_bootstrap_samples.R*: Bootstraps (and saves) m out of n ZIP codes.
  
  - *boot_associational.R*, *boot_gps_weighting_by_zip_year.R*, and *boot_gps_matching_by_zip_year.R*: R scripts for analyses on bootstrap samples.
  
  - *boot_associational.sh*, *boot_gps_weighting_by_zip_year.sh*, and *boot_gps_matching_by_zip_year.sh*: Bash script to run the aforementioned R scripts as a batch array job.
  
  - *process_parametric_boostrap_samples.R*: Collects the regression coefficient from each bootstrap analysis to calculate the m-out-of-n standard error.
  
`code/plot_results/`

  - *plot_covariate_balance_results.R*: Plots absolute correlation between exposure and covariates, for all exposures and all models.
    
  - *plot_parametric_results.R*: Plots hazard ratios for an IQR increase in each exposure, with 95% confidence intervals (m-out-of-n bootstrapped), for all exposures and all models.
  
  - *plot_smooths.R*: Plots natural cubic spline basis with 4 degrees of freedom, as a sensitivity analysis (alternative outcome model).
  
`code/exploratory/`

  - *explore_indiv_vars.R*: Creates Table 1.
    
  - *explore_zip_year_exposures.R*: Creates table of exposure summary statistics for the Appendix.

`results/`: Outputs saved by the above files.


**Package Versions**: All the analyses are run on R version 4.2.0 (2022-04-22) and [CausalGPS](https://github.com/cran/CausalGPS) version 0.3.0.9000.
