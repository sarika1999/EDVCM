# April 2024
# Run sampler 

# arguments from Rscript
args <- commandArgs(trailingOnly=TRUE)
sim_val = as.numeric(args[1]) #1-nsim

setwd('/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/') #cluster 

suppressMessages(library(tidyverse))

data <- readRDS('/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/data/general_simulation3/allfloodzips_dur3_mult_gp_general_simulation3.rds')
Y <- readRDS('/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/data/general_simulation3/allfloodzips_dur3_Y_general_simulation3.rds')

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
                     Sigma = data$Sigma,
                     sequence = data$sequence,
                     mu = data$mu)


library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = FALSE)

suppressMessages(
test <- stan(
  file = 'code/mult_gp_general_simulation3.stan',  # Stan program
  data = mult_gp_data,    # named list of data
  chains = 1, 
  iter = 500
))

#saveRDS(test, paste0('/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/output/simulations/general_simulation3/allfloodzips_dur3', sim_val, '.rds'))

beta_posterior <- as.data.frame(test) %>% select(-c("lp__"))

write.csv(beta_posterior, paste0('/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/output/simulations/general_simulation3/allfloodzips_dur3_', sim_val, '.csv'))
