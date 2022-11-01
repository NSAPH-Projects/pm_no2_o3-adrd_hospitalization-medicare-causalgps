Project: Causal_ADRD
Authors: Daniel Mork, Michelle Qin
Created: Jan. 21, 2022
Overview: Causal inference analysis of air pollution (PM2.5, NO2, Ozone)
  effect on AD/ADRD first hospitalization.
Data: all Medicare beneficiaries continuously enrolled in FFS plan, until
  first AD/ADRD hospitalization or change to HMO plan or death, aggregated
  by entry year, event year, zip, age, race, sex, dual eligible

File descriptions:
/R/
  /code/
    /aggregation/
      These five files create the dataset to be analyzed.

      1. Medicare FFS enrollment.R: reads denominator files, gets 2000-2016 FFS data, defines variables like race and sex as at entry into FFS
      2. ADRD hospitalization data.R: get first AD/ADRD hospitalization year, age, and zip code; merge into qid entry exit data
      3. Aggregate data sets.R: aggregate individuals by ffs_entry_year, AD_zip, AD_age_binned, AD_year, sexM, race_cat, any_dual
      4. Year-zip confounders.R: get confounders and exposures by zip and year
      5. Merge and clean: merge patient data with confounders and exposures, use complete cases
    
    /analysis/
      /exploratory/explore_covariate_distributions.R: examine distributions of covariates

      helper_functions.R: hard-code exposure name and covariates (classified as quantitative or categorical) as vectors, define some useful functions for transforming variables, checking covariate balance, etc.
      gps_models.R: model the GPS using CausalGPS package and other methods
      outcome_models.R: model the outcome using Poisson regression and (to do) semiparametric and nonparametric models

    

The following is not updated. For now, see files in code/aggregation/ to see what data files there are.

Data descriptions:
../../Data/Causal_ADRD/
  /aggregated/
    qid_agg_entry_to_followup.fst: complete merged data aggregated by zip,
      entry year, followup year, sex, race, dual enrollment, number AD/ADRD
      events, number deaths

  /denom/
    qid_complete.fst: individual data (1 row per unique individual-year, with
      full year enrollment) and individual specific data: zip, sex, age, race, 
      dual, death
    zip_yr_confounders.fst: confounders from census, bfrss, etc. at every zip 
      and year combination, ready for merging with aggregated obsevational data

