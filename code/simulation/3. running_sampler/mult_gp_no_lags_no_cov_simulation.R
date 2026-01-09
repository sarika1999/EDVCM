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
                     N = data$N,
                     rows_per_strata = data$rows_per_strata,
                     Y = Y[,sim_val],
                     X = data$X,
                     offset = data$offset,
                     Sigma_d = data$Sigma_d,
                     Sigma_t = data$Sigma_t,
                     sequence = data$sequence,
                     mu = data$mu)


library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = FALSE)

suppressMessages(
test <- stan(
  file = 'mult_gp_no_lags_no_cov_simulation.stan',  # Stan program
  data = mult_gp_data,    # named list of data
  chains = 4, 
  iter = 2000
))

saveRDS(test)

posterior <- as.data.frame(test) %>% select(-c("lp__")) 
names_col_beta <- paste0('beta',1:data$num_coeff)
colnames(posterior) <- c("little_sigma2", "phi", "tau", names_col_beta)

output_dir <- paste0('/n/holylabs/nethery_lab/Lab/floods-hospitalizations-glm/multinomial_GP/output/simulations/general_simulation', simulation, '/')
ifelse(!dir.exists(output_dir), dir.create(output_dir, recursive=TRUE), FALSE)

write.csv(posterior, 
          paste0(output_dir, 'rand_general_smooth_med_', sim_val, '.csv'),
          row.names=FALSE)


