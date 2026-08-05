# June 2024
# Run sampler 

# arguments from Rscript
#args <- commandArgs(trailingOnly=TRUE)

setwd('./multinomial_GP/') #cluster 

suppressMessages(library(tidyverse))

data_dir <- './multinomial_GP/FFS_final_oct2025/sampler/'

#note this file does not exist currently and would need to be created#
data <- readRDS(paste0(data_dir, 'cause_2000_2016_spline_covariates_lags_dur10_mult_gp_real_data.rds'))

mult_gp_data <- list(floodzip_id = data$floodzip_id,
                     case_control_set = data$case_control_set,
                     durations = data$durations,
                     D = data$D,
                     num_coeff = data$num_coeff,
                     lag = data$lag,
                     num_lag_coeff = data$num_lag_coeff,
                     N = data$N,
                     rows_per_strata = data$rows_per_strata,
                     Y = Y[,sim_val],
                     A = data$A,
                     X = data$X,
                     offset = data$offset,
                     Sigma_d = data$Sigma_d,
                     Sigma_t = data$Sigma_t,
                     Sigma_d2 = data$Sigma_d2,
                     Sigma_l = data$Sigma_l,
                     sequence = data$sequence,
                     mu_beta = data$mu_beta,
                     mu_theta = data$mu_theta)

library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

suppressMessages(
  test <- stan(
    file = 'code/EDVCM_cov_lags_application.stan',  # Stan program
    data = mult_gp_data,    # named list of data
    chains = 1,
    iter = 2000,
    warmup = 500,
  ))

output_dir <- './multinomial_GP/FFS_final_oct2025/results/'
ifelse(!dir.exists(output_dir), dir.create(output_dir, recursive=TRUE), FALSE)

saveRDS(test, paste0(output_dir, 'fit_cause_with_lags.rds'))

posterior <- as.data.frame(test) %>% select(-c("lp__")) 
names_col_beta <- paste0('beta',1:data$num_coeff)
names_col_theta <- paste0('theta', 1:data$num_lag_coeff)
names_col_conf <- paste0('conf', 1:(data$num_conf*data$spline_df))
colnames(posterior) <- c("sigma_beta", "phi", "tau", "sigma_theta", "gamma", "eta", names_col_beta, names_col_theta, names_col_conf, "little_sigma_beta2", "little_sigma_theta2")


write.csv(posterior,
          paste0(output_dir, 'cause_2000_2016_covariates_lags_dur10.csv'),
          row.names=FALSE)

