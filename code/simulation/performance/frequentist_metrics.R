# frequentist performance metrics

setwd('./multinomial_GP/') 

# Load necessary library
library(dplyr)
library(purrr)

# Set the directory containing the CSV files
data_dir <- './multinomial_GP/output/simulations/freq/med_smooth_spline_Xpctnoise/results/unconstrained_dlm/summary/'

combined_data <- vector("list", length = 5000) #5000 is n_sim
for (i in 1:5000){
  # Get a list of all CSV files in the directory
  csv_files <- list.files(path = data_dir, pattern = paste0("rand_data", i,"_.*\\.csv"), full.names = TRUE)
  
  # Read and combine all CSV files into one data frame
  combined_data[[i]] <- csv_files %>%
    lapply(read.csv) %>%
    bind_rows()
  
  combined_data_cleaned <- lapply(combined_data, function(df) df[c(1, 62:105, 2:61), ]) #reorder to be in duration order #depends on D
  
  # View the combined data
  head(combined_data_cleaned)
}

saveRDS(combined_data_cleaned, "freq_Xpctnoise_combined_results.rds")

combined_data_cleaned_freqX <- readRDS('./simulations/freq/med_smooth_spline_Xpctnoise/freq_Xpctnoise_combined_results.rds')

# combined_data_cleaned: list of 5000 data.frames
# each data.frame has columns: X, t, rr, rr.ll, rr.ul

true_beta <- readRDS('./multinomial_GP/data/simulations/general_simulation#/rand_true_beta_general_smooth_med_Xpctnoise_simulation#.rds')

# Extract third row from each data frame
third_rows <- lapply(combined_data_cleaned, function(df) df[,3])

# Combine extracted rows into a matrix
result_matrix <- do.call(rbind, third_rows)

# Mean estimate per coefficient across simulations
mean_estimates <- colMeans(result_matrix)

# Compute percent bias
percent_bias <- ((mean_estimates - true_beta) / true_beta) * 100

saveRDS(percent_bias, 'percent_bias_freqX.rds') #latest version is saved as JUN25

# Create a matrix of true values repeated across simulations
repeat_mat <- matrix(true_beta, nrow = nrow(result_matrix), ncol = length(true_beta), byrow = TRUE)

# Compute MSE per coefficient
mse <- colMeans((result_matrix - repeat_mat)^2)

saveRDS(mse, 'mse_freqX.rds')

#coverage 

# For each simulation, compute a logical vector:
#    TRUE if rr is inside [rr.ll, rr.ul], FALSE otherwise
coverage_mat_freqX <- sapply(combined_data_cleaned_freqX, function(d) {
  (rand_true_beta_general_smooth_med_Xpctnoise_simulationY >= d$rr.ll & rand_true_beta_general_smooth_med_Xpctnoise_simulationY <= d$rr.ul)
})

# coverage_mat will be a ((D*D+1)/2) x n_sim logical matrix
dim(coverage_mat_freqX) 

# Coverage probability for each row (parameter)
coverage_prob_freqX <- rowMeans(coverage_mat_freqX)
saveRDS(coverage_prob_freqX, 'coverage_prob_freqX.rds')

