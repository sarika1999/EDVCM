# April 2024
# Check bias and coverage 

simulation <- 4
output_dir <- setwd(paste0('/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/output/simulations/general_simulation', simulation, '/')) #cluster 
data_dir <- paste0('/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/data/simulations/general_simulation', simulation, '/')

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

true_beta <- readRDS(paste0(data_dir, 'allfloodzips_dur3_true_beta_general_simulation', simulation, '.rds'))
hyperparameter <- c(0.2, 0.5, 0.5)

posterior_list <- list.files(pattern="\\.csv$")
posterior <- lapply(posterior_list, fread)

#bias
posterior_bias <- lapply(posterior, mutate, sigma_bias = little_sigma2 - hyperparameter[1],
                                                      phi_bias = phi - hyperparameter[2],
                                                      tau_bias = tau - hyperparameter[3],
                                                      bias1 = beta1 - true_beta[1],
                                                      bias2 = beta2 - true_beta[2],
                                                      bias3 = beta3 - true_beta[3],
                                                      bias4 = beta4 - true_beta[4],
                                                      bias5 = beta5 - true_beta[5],
                                                      bias6 = beta6 - true_beta[6])

posterior_mean_bias <- matrix(data = NA, nrow = length(posterior_bias), ncol = (length(hyperparameter) + length(true_beta)))
for (sim in 1:length(posterior_bias)){
  posterior_mean_bias[sim,] <- colMeans(posterior_bias[[sim]])[7:15]
}

boxplot(posterior_mean_bias)

#coverage 
posterior_quantiles <- matrix(data = NA, nrow = length(posterior), ncol = (2*(length(hyperparameter) + length(true_beta))))
for (sim in 1:length(posterior)){
  posterior_quantiles[sim,] <- apply(posterior[[sim]], 2, quantile, probs = c(0.025, 0.975))
}

colnames(posterior_quantiles) <- c("littlesigma2_l", "littlesigma2_u", "phi_l","phi_u", "tau_l", "tau_u", 
                                   "beta1_l", "beta1_u", "beta2_l", "beta2_u", "beta3_l", "beta3_u",
                                   "beta4_l", "beta4_u", "beta5_l", "beta5_u", "beta6_l", "beta6_u")

coverage <- matrix(data = NA, nrow = length(posterior), ncol = (length(hyperparameter) + length(true_beta)))

coverage[,1] <- ifelse(posterior_quantiles[,c(1)] <= hyperparameter[1] & posterior_quantiles[,c(2)] >= hyperparameter[1], 1, 0)
coverage[,2] <- ifelse(posterior_quantiles[,c(3)] <= hyperparameter[2] & posterior_quantiles[,c(4)] >= hyperparameter[2], 1, 0)
coverage[,3] <- ifelse(posterior_quantiles[,c(5)] <= hyperparameter[3] & posterior_quantiles[,c(6)] >= hyperparameter[3], 1, 0)
coverage[,4] <- ifelse(posterior_quantiles[,c(7)] <= true_beta[1] & posterior_quantiles[,c(8)] >= true_beta[1], 1, 0)
coverage[,5] <- ifelse(posterior_quantiles[,c(9)] <= true_beta[2] & posterior_quantiles[,c(10)] >= true_beta[2], 1, 0)
coverage[,6] <- ifelse(posterior_quantiles[,c(11)] <= true_beta[3] & posterior_quantiles[,c(12)] >= true_beta[3], 1, 0)
coverage[,7] <- ifelse(posterior_quantiles[,c(13)] <= true_beta[4] & posterior_quantiles[,c(14)] >= true_beta[4], 1, 0)
coverage[,8] <- ifelse(posterior_quantiles[,c(15)] <= true_beta[5] & posterior_quantiles[,c(16)] >= true_beta[5], 1, 0)
coverage[,9] <- ifelse(posterior_quantiles[,c(17)] <= true_beta[6] & posterior_quantiles[,c(18)] >= true_beta[6], 1, 0)

coverage_prob <- colMeans(coverage)
coverage_prob






