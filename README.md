# A varying-coefficient model for characterizing duration-driven heterogeneity in flood-related health impacts

This repository contains R and corresponding Stan code for the simulations and applications presented in manuscript: 

## Project Structure

### Data 
#### `simulation/`: This directory details simulation data that was used to run the simulation study.
FASRC location: ./multinomial_GP/
- 45: medium smoothing spline + 25\% random noise (GP prior) 
- 48: medium smoothing spline + 100\% random noise (GP prior) 
- 49: medium smoothing spline + 25\% random noise (same underlying surface values as 45 but ran with N(0,1) prior) 
- 52: medium smoothing spline + 100\% random noise (same underlying surface values as 48 but ran with N(0,1) prior)
- 53: medium smoothing spline + 25\% random noise + missing durations (same underlying surface vales as 45 but removed durations, ran with GP prior) 
- 54: medium smoothing spline + 100\% random noise + missing durations (same underlying surface vales as 45 but removed durations, ran with GP prior) 
- 56: medium smoothing spline + 25\% random noise for beta (same underlying surface values as 48), medium smoothing spline + 25\% random noise for theta 
- 57: medium smoothing spline + 100\% random noise for beta (same underlying surface values as 48), medium smoothing spline + 100\% random noise for theta 
- frequent_comparison_Xpctnoise: distinct datasets for each duration for each simulation iterate with the same outcome counts as 45, 48 

Notes
1. Medium smoothing spline refers to a thin-plate spline with 5 basis functions. 
2. Simulations rely on the following data structures: 
    - rand_data.rds was created from events_with_matched_controls_nolag_df_v3.rds
    - 45, 49, 48, 52: rand_data.rds, X_randdata.rds
    - 53, 54: rand_data_no_dur4_dur7_dur11.rds, X_randdata_no_dur4_dur7_dur11.rds
    - 56, 57: lags_rand_data.rds, X_lags_2d_randdata.rds, A_lags_2d_randdata.rds 
3. In 56 and 57, the underlying smoothing spline for the lags is the same (true_theta.rds) and then varying amounts of random noise is added (similar to the primary effect construction). 
4. Simulations include durations 1-14 and 5 (2-dimensional) lagged effects for each duration when applicable. 

#### `application/`: This directory details data used for an application on floods and cause-specific hospitalization.
ReD location: ./multinomial_GP/FFS_final_oct2025

### Code 

#### `simulation/`: This directory includes all code used to run the simulation study.
1. dataset creation 
- `random_data_sample.R`
- `make_lags_data.R`
2. setup 
- `smoothing_simulation_setup.R`
- `smoothing_simulation_setup_lags.R`
- `ground_truth_heatmap.R`
3. running sampler 
- `mult_gp_no_lags_no_cov_simulation.R`
- `mult_gp_no_lags_no_cov_simulation.stan`
- `mult_gp_no_lags_no_cov_simulation_comparator.R`
- `mult_gp_no_lags_no_cov_simulation_comparator.stan`
- `mult_gp_lags_2d_simulation.R`
- `mult_gp_lags_2d_simulation.stan`

#### `application/`: This directory includes example code for running the real data application. 

#### `performance/`: This directory consists of code used to evaluate model performance from the simulation study including figure generation. 
- `metrics.R`: compute percent bias, mean squared error, and coverage 
- `frequentist_metrics.R`: compute percent bias, mean squared error, and coverage for the frequentist comparison
- `additional_metrics.R`: obtain treedepth, effective sample size, and rhat; obtain significance using 95\% credible interval and direction for significant coefficient estimates
- `heatmaps.R`: plot metrics for each duration-day coefficient estimate (or lagged day) 


## Usage Instructions

## Dependencies

The project requires the following dependencies:
- R 4.2.0+
- rstan
- tidyverse
- ggplot2
- dplyr
- tidyr
- purrr
- magrittr
- knitr
- rmarkdown

