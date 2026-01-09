
simulation <- 
setwd(paste0('/n/holylabs/nethery_lab/Lab/floods-hospitalizations-glm/multinomial_GP/output/simulations/general_simulation', simulation, '/')) #cluster
output_dir <- paste0('/n/holylabs/nethery_lab/Lab/floods-hospitalizations-glm/multinomial_GP/output/simulations/general_simulation', simulation, '/')
#figures_dir <- paste0('/n/holylabs/nethery_lab/Lab/floods-hospitalizations-glm/multinomial_GP/output/simulations/general_simulation', simulation, '/figures/')
#ifelse(!dir.exists(figures_dir), dir.create(figures_dir, recursive=TRUE), FALSE)
data_dir <- paste0('/n/holylabs/nethery_lab/Lab/floods-hospitalizations-glm/multinomial_GP/data/simulations/general_simulation', simulation, '/')

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

true_beta <- readRDS(paste0(data_dir, 'rand_true_beta_general_smooth_med_.rds')) #primary effects 
true_theta <- readRDS(paste0(data_dir, 'rand_true_theta_general_smooth_med_.rds')) #lags
hyperparameter <- NULL #we do not set these a priori when we use the smoothing spline 

parameters <- c(hyperparameter, true_beta, true_theta)

#some replicates may not finish due to time
posterior_list <- list.files(pattern="\\.csv$")

posterior <- lapply(posterior_list, fread)
#start at 4 if no lags (3 hyperparameters) #start at 7 if 2d lags (6 hyperparameters) #start at 5 if 1d lags or if 2d lags with shared hyperparameter for duration
posterior <- lapply(posterior, function(x) { x[,7:ncol(x)] }) 

# posterior = list(matrix1, matrix2, ..., matrix)
# matrix: n_sim * length(true_beta) = 5000 x 105
# true_beta: length D*(D+1)/2 = 105 
# true_theta: length l (1d lags) or l*D (2d lags)

# Step 1: Compute posterior means for each simulation
posterior_means <- lapply(posterior, colMeans)  # list of 5000 vectors, each of length 105

# Step 2: Convert to a matrix: 5000 x 105
posterior_means_mat <- do.call(rbind, posterior_means)

# Step 3: Compute average bias and MSE for each coefficient
bias <- colMeans(posterior_means_mat) - c(true_beta, true_theta)

percent_bias <- (bias / c(true_beta, true_theta)) * 100

mse <- colMeans((posterior_means_mat - matrix(parameters, nrow = nrow(posterior_means_mat), ncol = ncol(posterior_means_mat), byrow = TRUE))^2)

saveRDS(percent_bias, "percent_bias.rds") #latest version is saved as JUN25

saveRDS(mse, "mse.rds")

#coverage
posterior_quantiles <- matrix(data = NA, nrow = length(posterior), ncol = (2*(length(parameters))))
for (sim in 1:length(posterior)){
  posterior_quantiles[sim,] <- apply(posterior[[sim]], 2, quantile, probs = c(0.025, 0.975))
}

colnames(posterior_quantiles) <- rep(colnames(posterior[[1]]), each = 2)
for (i in 1:ncol(posterior_quantiles)){
  if (i %% 2 == 0){
    colnames(posterior_quantiles)[i] <- paste0(colnames(posterior_quantiles)[i], "_u")
  }
  else {
    colnames(posterior_quantiles)[i] <- paste0(colnames(posterior_quantiles)[i], "_l")
  }
}

#get column position for lower and upper bound of each parameter
meder <- seq(1,2*length(parameters), 2)
upper <- seq(2,2*length(parameters), 2)

coverage <- matrix(data = NA, nrow = length(posterior), ncol = length(parameters))
for (i in 1:length(parameters)){
  coverage[,i] <- ifelse(posterior_quantiles[,meder[i]] <= parameters[i] & posterior_quantiles[,upper[i]] >= parameters[i], 1, 0)
}

coverage_prob <- colMeans(coverage)
saveRDS(coverage_prob, paste0(output_dir, "coverage_prob.rds"))

