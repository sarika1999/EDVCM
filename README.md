# A varying-coefficient model for characterizing duration-driven heterogeneity in flood-related health impacts

This repository contains R and corresponding Stan code for the simulations and applications presented in manuscript: 

## Project Structure

### Data 
#### `simulation/`: This directory details simulation data that was used to run the simulation study.
FASRC location: /n/holylabs/nethery_lab/Lab/floods-hospitalizations-glm/multinomial_GP/data/simulations
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

#### `simulation/`: This directory includes all code used to run and evaluate the simulation study including figure generation.
1. dataset creation 
- `random_data_sample.R`: take 10\% of real exposure dataset to use in simulation study 
- `make_lags_data.R`: add lagged days (post-flood days) to dataset 
2. setup 
- `smoothing_simulation_setup.R`: 
- `smoothing_simulation_setup_lags.R`:
- `ground_truth_heatmap.R`: make heatmap of true underlying surfaces (generated coefficient values) 
3. running sampler 
- `mult_gp_no_lags_no_cov_simulation.R`:
- `mult_gp_no_lags_no_cov_simulation.stan`:
- `mult_gp_no_lags_no_cov_simulation_comparator.R`:
- `mult_gp_no_lags_no_cov_simulation_comparator.stan`:
- `mult_gp_lags_2d_simulation.R`:
- `mult_gp_lags_2d_simulation.stan`:
- `frequentist_simulation.R`: 
4. performance
- `metrics.R`: compute percent bias, mean squared error, and coverage 
- `frequentist_metrics.R`: compute percent bias, mean squared error, and coverage for the frequentist comparison
- `additional_metrics.R`: obtain treedepth, effective sample size, and rhat
- `metric_heatmaps.R`: make heatmap of each metric for each duration-day (or lagged day) coefficient point estimate 

#### `application/`: This directory includes example code for running the real data application including figure generation. 
1. dataset creation 
- `make_data.R`: 
- `combine_county_data_cov_lags.R`: 
2. setup 
- `cov_spline_setup.R`: create spline basis matrix for each covariate to input into stan sampler 
- `application_setup.R`: take all data pieces and output a single list to input into stan sampler 
3. running sampler 
- `mult_gp_cov_application.R`: R script for running sampler on a particular cause of hospitalization
- `mult_gp_cov_application.stan`: stan sampler for application (includes time-varying covariates)
4. results
- `model_fit.R`: obtain point estimates and corresponding 95\% credible intervals; determine significance and direction for significant coefficient estimates; format point estimates and credible intervals for covariate terms
- `cumulative_effects.R`: calculate cumulative rate ratio for each duration with(out) presence of time-varying covariates 
- `point_estimate_heatmap.R`: make heatmap of duration-day coefficient point estimates (rate ratios) 

#### `misc/`: This directory includes miscellaneous code used to create figures that are not specific to the method.
- `spatial_trends.R`: make county-level map of frequency (number of times flooded) and average duration over the study period 

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

