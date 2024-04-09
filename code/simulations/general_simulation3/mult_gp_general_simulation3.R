# April 2024
# Run sampler 

# arguments from Rscript
args <- commandArgs(trailingOnly=TRUE)
sim_val = as.numeric(args[1]) #1-nsim

setwd('/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/') #cluster 

suppressMessages(library(tidyverse))

data <- readRDS('/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/data/general_simulation3/allfloodzips_dur3_mult_gp_general_simulation3.rds')
Y <- readRDS('/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/data/general_simulation3/allfloodzips_dur3_Y_general_simulation3.rds')

mult_gp_data <- list(floodzip_id = floodzip_id,
                     case_control_set = case_control_set,
                     durations = durations,
                     D = D,
                     num_coeff = num_coeff,
                     N = N,
                     rows_per_strata = rows_per_strata,
                     Y = Y[,sim_val],
                     X = X,
                     offset = offset,
                     Sigma = Sigma,
                     sequence = sequence,
                     mu = mu)


library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

suppressMessages(
  test <- stan(
    file = 'code/mult_gp_general_simulation3.stan',  # Stan program
    data = mult_gp_data,    # named list of data
    chains = 4, 
    iter = 2000
))

#saveRDS(test, paste0('/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/output/simulations/general_simulation3/allfloodzips_dur3', sim_val, '.rds'))

beta_posterior <- as.data.frame(test) %>% select(-c("lp__")) %>% rename("beta[1][1]" = "beta[1]",
                                                                        "beta[2][1]" = "beta[2]",
                                                                        "beta[2][2]" = "beta[3]",
                                                                        "beta[3][1]" = "beta[4]",
                                                                        "beta[3][2]" = "beta[5]",
                                                                        "beta[3][3]" = "beta[6]")

write.csv(beta_posterior, paste0('/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/output/simulations/general_simulation3/allfloodzips_dur3', sim_val, '.csv'))

