# April 2024
# Run sampler 

# arguments from Rscript
args <- commandArgs(trailingOnly=TRUE)
sim_val = as.numeric(args[1]) #1-nsim

setwd('/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/') #cluster 

suppressMessages(library(tidyverse))

simulation <- 13
data_dir <- paste0('/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/data/simulations/general_simulation', simulation, '/')

data <- readRDS(paste0(data_dir, 'rand_lnormprior_mult_gp_general_simulation', simulation, '.rds'))
Y <- readRDS(paste0(data_dir, 'rand_lnormprior_Y_general_simulation', simulation, '.rds'))

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
  file = 'code/mult_gp_general_simulation_dupe.stan',  # Stan program
  data = mult_gp_data,    # named list of data
  chains = 2, 
  iter = 2000
))

#saveRDS(test, paste0('/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/output/simulations/general_simulation3/allfloodzips_no_dur1_', sim_val, '.rds'))

posterior <- as.data.frame(test) %>% select(-c("lp__")) 
#colnames(posterior) <- c("little_sigma2", "phi", "tau", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6")
names_col_beta <- paste0('beta',1:data$num_coeff)
colnames(posterior) <- c("little_sigma2", "phi", "tau", names_col_beta)

output_dir <- paste0('/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/output/simulations/general_simulation', simulation, '/')
ifelse(!dir.exists(output_dir), dir.create(output_dir, recursive=TRUE), FALSE)

write.csv(posterior, 
          paste0(output_dir, 'rand_lnormprior_', sim_val, '.csv'),
          row.names=FALSE)


