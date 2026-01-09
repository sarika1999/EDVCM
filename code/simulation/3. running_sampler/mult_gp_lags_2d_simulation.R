# Run sampler 

# arguments from Rscript
args <- commandArgs(trailingOnly=TRUE)
sim_val = as.numeric(args[1]) #1-nsim

setwd('/n/holylabs/nethery_lab/Lab/floods-hospitalizations-glm/multinomial_GP/') #cluster 

suppressMessages(library(tidyverse))

simulation <- 
data_dir <- paste0('/n/holylabs/nethery_lab/Lab/floods-hospitalizations-glm/multinomial_GP/data/simulations/general_simulation', simulation, '/')

data <- readRDS(paste0(data_dir, 'rand_mult_gp_general_smooth_med_', simulation, '.rds'))
Y <- readRDS(paste0(data_dir, 'rand_Y_general_smooth_med_', simulation, '.rds'))

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
rstan_options(auto_write = FALSE)

suppressMessages(
test <- stan(
  file = 'mult_gp_lags_2d_simulation.stan',  # Stan program
  data = mult_gp_data,    # named list of data
  chains = 4, 
  iter = 2000
))

saveRDS(test)

posterior <- as.data.frame(test) %>% select(-c("lp__")) 
names_col_beta <- paste0('beta',1:data$num_coeff)
names_col_theta <- paste0('theta', 1:data$num_lag_coeff)
colnames(posterior) <- c("little_sigma_beta2", "phi", "tau", "little_sigma_theta2", "gamma", "eta", names_col_beta, names_col_theta)

output_dir <- paste0('/n/holylabs/nethery_lab/Lab/floods-hospitalizations-glm/multinomial_GP/output/simulations/general_simulation', simulation, '/')
ifelse(!dir.exists(output_dir), dir.create(output_dir, recursive=TRUE), FALSE)

write.csv(posterior, 
          paste0(output_dir, 'rand_lags_2d_', sim_val, '.csv'),
          row.names=FALSE)


