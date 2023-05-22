**Project**: Causal_ADRD

**Authors**: Daniel Mork, Michelle Qin, Naeem Khoshnevis

**Created**: Jan. 21, 2022

**Overview**: Causal inference analysis of long-term air pollution (PM2.5, NO2, summer ozone) effect on first hospitalization with Alzheimer's disease and related dementias (ADRD).

**Data**: All Medicare beneficiaries continuously enrolled in FFS plan between 2001 and 2016, until first ADRD hospitalization or change to HMO plan or death, aggregated by entry year, event year, ZIP code, 5-year-binned age group, race, sex, and Medicaid eligibility.

**File and folder descriptions**: 

`aggregation`   

These six files create the three datasets, one for each exposure, to be analyzed.   

- *1. Medicare FFS enrollment.R*: Reads denominator files, gets 2000-2016 FFS data, defines variables like race and sex as at entry into FFS.
- *2. ADRD hospitalization data.R*: Gets first ADRD hospitalization year, age, and ZIP code; merge into qid entry exit data.
- *3. Aggregate data sets.R*: Aggregates individuals by FFS entry year (called `cohort`), `zip`, `age_grp`, `year`, `sex`, `race`, and Medicaid eligibility (called `dual`). Each observation is now an "aggregated patient unit".
- *4. Year-zip confounders.R:* Gets confounders and exposures by zip and year.
- *5. Merge and clean.R*: Merges patient data with confounders from current year and exposures from previous year, using complete cases.
- *6. Set up trimmed and ZIP data.R*: For each of PM2.5, NO2, summer ozone at the ZIP code and annual level, trims the dataset at the exposure's 5th and 95th percentiles, keeping the middle 90\% of aggregated patient units.
    
`analysis`

- *exploratory/explore_covariate_distributions.R*: examine distributions of covariates    
- *helper_functions.R*: hard-code exposure name and covariates (classified as quantitative or categorical) as vectors, define some useful functions for transforming variables, checking covariate balance, etc.   
- *gps_models.R*: model the GPS using CausalGPS package and other methods
- *outcome_models.R*: model the outcome using Poisson regression and (to do) semiparametric and nonparametric models
