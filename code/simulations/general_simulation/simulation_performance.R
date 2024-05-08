# April 2024
# Check bias and coverage 

simulation <- 9
output_dir <- setwd(paste0('/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/output/simulations/general_simulation', simulation, '/')) #cluster 
figures_dir <- paste0('/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/output/simulations/general_simulation', simulation, '/figures/')
ifelse(!dir.exists(figures_dir), dir.create(figures_dir, recursive=TRUE), FALSE)
data_dir <- paste0('/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/data/simulations/general_simulation', simulation, '/')

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

true_beta <- readRDS(paste0(data_dir, 'rand_no_dur2_true_beta_general_simulation', simulation, '.rds'))
hyperparameter <- c(0.2, 0.4, 0.6)

parameters <- c(hyperparameter, true_beta)

#some replicates may not finish due to time 
posterior_list <- list.files(pattern="\\.csv$")
posterior <- lapply(posterior_list, fread)

# bias
compute_bias <- function(df, truth) {
  #for each column (parameter) in the data frame, compute the |sample - true| 
  df <- map_dfc(1:ncol(df), function(i) {
    df_tmp <- tibble('tmp' = abs(df[[i]] - truth[i]))
    #make the colname of the new data frame the same as before (parameters)
    names(df_tmp) <- names(df)[i]
    return(df_tmp)
  })
  return(df)
}

library(furrr) #allows 'map' to run in parallel 
posterior_bias <- future_map(posterior, ~compute_bias(.x, parameters)) #apply function above to all data frames in the list 

#compute the mean bias for each parameter (for each data frame in list)
posterior_mean_bias <- matrix(data = NA, nrow = length(posterior_bias), ncol = length(parameters))
for (sim in 1:length(posterior_bias)){
  posterior_mean_bias[sim,] <- colMeans(posterior_bias[[sim]])
}

boxplot(posterior_mean_bias)

summary_mean_bias <- matrix(data = NA, nrow = 2, ncol = ncol(posterior_mean_bias))
for (i in 1:ncol(posterior_mean_bias)){
  summary_mean_bias[1,col] <- median(posterior_mean_bias[,i])
  summary_mean_bias[2,col] <- max(posterior_mean_bias[,i])
}

which.max(summary_mean_bias[1,4:108])

#coverage 
posterior_quantiles <- matrix(data = NA, nrow = length(posterior), ncol = (2*length(parameters)))
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
lower <- seq(1,2*length(parameters), 2)
upper <- seq(2,2*length(parameters), 2)

coverage <- matrix(data = NA, nrow = length(posterior), ncol = length(parameters))
for (i in 1:length(parameters)){
  coverage[,i] <- ifelse(posterior_quantiles[,lower[i]] <= parameters[i] & posterior_quantiles[,upper[i]] >= parameters[i], 1, 0)
}

coverage_prob <- colMeans(coverage)





