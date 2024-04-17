# April 2024
# Check bias and coverage 

setwd('/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/output/simulations/general_simulation3/') #cluster 

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

true_beta <- readRDS('/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/data/general_simulation3/allfloodzips_dur3_true_beta_general_simulation3.rds')

beta_list <- list.files(pattern="\\.csv$")
beta_posterior <- lapply(beta_list, fread)

#bias
beta_posterior_bias <- lapply(beta_posterior, mutate, bias1 = beta1 - true_beta[1],
                                                         bias2 = beta2 - true_beta[2],
                                                         bias3 = beta3 - true_beta[3],
                                                         bias4 = beta4 - true_beta[4],
                                                         bias5 = beta5 - true_beta[5],
                                                         bias6 = beta6 - true_beta[6])

beta_posterior_mean_bias <- matrix(data = NA, nrow = length(beta_posterior_bias), ncol = length(true_beta))
for (sim in 1:length(beta_posterior_bias)){
  beta_posterior_mean_bias[sim,] <- colMeans(beta_posterior_bias[[sim]])[7:12]
}

boxplot(beta_posterior_mean_bias)

#coverage 
beta_posterior_quantiles <- matrix(data = NA, nrow = length(beta_posterior), ncol = 2*length(true_beta))
for (sim in 1:length(beta_posterior)){
  beta_posterior_quantiles[sim,] <- apply(beta_posterior[[sim]], 2, quantile, probs = c(0.025, 0.975))
}

colnames(beta_posterior_quantiles) <- c("beta1_l", "beta1_u", "beta2_l", "beta2_u", "beta3_l", "beta3_u",
                                        "beta4_l", "beta4_u", "beta5_l", "beta5_u", "beta6_l", "beta6_u")

coverage <- matrix(data = NA, nrow = length(beta_posterior), ncol = length(true_beta))

coverage[,1] <- ifelse(beta_posterior_quantiles[,c(1)] <= true_beta[1] & beta_posterior_quantiles[,c(2)] >= true_beta[1], 1, 0)
coverage[,2] <- ifelse(beta_posterior_quantiles[,c(3)] <= true_beta[2] & beta_posterior_quantiles[,c(4)] >= true_beta[2], 1, 0)
coverage[,3] <- ifelse(beta_posterior_quantiles[,c(5)] <= true_beta[3] & beta_posterior_quantiles[,c(6)] >= true_beta[3], 1, 0)
coverage[,4] <- ifelse(beta_posterior_quantiles[,c(7)] <= true_beta[4] & beta_posterior_quantiles[,c(8)] >= true_beta[4], 1, 0)
coverage[,5] <- ifelse(beta_posterior_quantiles[,c(9)] <= true_beta[5] & beta_posterior_quantiles[,c(10)] >= true_beta[5], 1, 0)
coverage[,6] <- ifelse(beta_posterior_quantiles[,c(11)] <= true_beta[6] & beta_posterior_quantiles[,c(12)] >= true_beta[6], 1, 0)

coverage_prob <- colMeans(coverage)
coverage_prob






