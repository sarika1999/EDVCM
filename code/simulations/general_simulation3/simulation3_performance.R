# April 2024
# Check bias and coverage 

setwd('/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/') #cluster 

suppressMessages(library(tidyverse))

true_beta <- readRDS('/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/data/general_simulation3/allfloodzips_dur3_true_beta_general_simulation3.rds')

for (sim_val in 1:1000){
  beta_posterior_list <- read.csv()
}

