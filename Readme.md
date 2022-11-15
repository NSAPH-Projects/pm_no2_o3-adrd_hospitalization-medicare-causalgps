**Project**: Causal_ADRD    
**Authors**: Daniel Mork, Michelle Qin    
**Created**: Jan. 21, 2022    
**Overview**: Causal inference analysis of air pollution (PM2.5, NO2, Ozone)
  effect on AD/ADRD first hospitalization.    
**Data**: all Medicare beneficiaries continuously enrolled in FFS plan, until
  first AD/ADRD hospitalization or change to HMO plan or death, aggregated
  by entry year, event year, zip, age, race, sex, dual eligible    

**File and folder descriptions**: 

`aggregation`   

These five files create the dataset to be analyzed.   

1. Medicare FFS enrollment.R: reads denominator files, gets 2000-2016 FFS data, defines variables like race and sex as at entry into FFS   
2.ADRD hospitalization data.R: get first AD/ADRD hospitalization year, age, and zip code; merge into qid entry exit data.  
3. Aggregate data sets.R: aggregate individuals by ffs_entry_year, AD_zip, AD_age_binned, AD_year, sexM, race_cat, any_dual.  
4. Year-zip confounders.R: get confounders and exposures by zip and year  
5. Merge and clean: merge patient data with confounders and exposures, use complete cases.  
    
`analysis`
- */exploratory/explore_covariate_distributions.R*: examine distributions of covariates    
- *helper_functions.R*: hard-code exposure name and covariates (classified as quantitative or categorical) as vectors, define some useful functions for transforming variables, checking covariate balance, etc.   
- *gps_models.R*: model the GPS using CausalGPS package and other methods
- *outcome_models.R*: model the outcome using Poisson regression and (to do) semiparametric and nonparametric models
