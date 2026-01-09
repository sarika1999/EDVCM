# June 2024
# Run sampler 

# arguments from Rscript
#args <- commandArgs(trailingOnly=TRUE)

setwd('./multinomial_GP/') #cluster 

suppressMessages(library(tidyverse))

data_dir <- './multinomial_GP/FFS_final_oct2025/sampler/'

data <- readRDS(paste0(data_dir, 'cause_2000_2016_spline_covariates_no_lag_dur10_mult_gp_real_data.rds'))

mult_gp_data <- list(floodcty_id = data$floodcty_id,
                     case_control_set = data$case_control_set,
                     durations = data$durations,
                     D = data$D,
                     num_coeff = data$num_coeff,
                     num_conf = data$num_conf, 
                     spline_df = data$spline_df,
                     N = data$N,
                     rows_per_strata = data$rows_per_strata,
                     Y = data$Y,
                     X = data$X,
                     Z = data$Z,
                     offset = data$offset,
                     Sigma_d = data$Sigma_d,
                     Sigma_t = data$Sigma_t,
                     sequence = data$sequence,
                     mu_beta = data$mu_beta)


library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

suppressMessages(
  test <- stan(
    file = 'code/mult_gp_cov_application.stan',  # Stan program
    data = mult_gp_data,    # named list of data
    chains = 1,
    iter = 2000,
    warmup = 500,
  ))

output_dir <- './multinomial_GP/FFS_final_oct2025/results/'
ifelse(!dir.exists(output_dir), dir.create(output_dir, recursive=TRUE), FALSE)

saveRDS(test, paste0(output_dir, 'fit_cause.rds'))

posterior <- as.data.frame(test) %>% select(-c("lp__")) 
names_col_beta <- paste0('beta',1:data$num_coeff)
#names_col_theta <- paste0('theta', 1:data$lag)
names_col_conf <- paste0('conf', 1:(data$num_conf*data$spline_df))
colnames(posterior) <- c("little_sigma_beta2", "phi", "tau", names_col_beta, names_col_conf) #sigma, phi, tau, beta, zeta, sigma^2


write.csv(posterior, 
          paste0(output_dir, 'cause_2000_2016_covariates_dur10.csv'),
          row.names=FALSE)

