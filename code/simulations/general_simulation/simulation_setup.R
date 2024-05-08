# April 2024
# Generate data from the true DGP (one realization of the GP prior)

setwd('/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/') #cluster 

suppressMessages(library(tidyverse))

data_setup <- function(path_to_data, little_sigma2, phi, tau){
  dat <- readRDS(paste0(path_to_data, '/test_data.rds'))
  dat <- dat %>% arrange(floodzip_id,zipcode,month,day,control_indicator) #order to match paper 
  floodzip_id <- length(unique(dat$floodzip_id)) %>% as.double()
  case_control_set <- length(unique(dat$control_indicator)) %>% as.double()
  durations <- dat %>% group_by(floodzip_id) %>% summarise(durations = unique(d)) %>% select(durations) %>% unname() %>% unlist()
  D <- max(durations)
  num_coeff <- D*(D+1)/2
  
  #define N as 0, then compute since there are multiple sets of rows per flood-zip combination if duration > 1
  #N <- nrow(dat)
  N <- 0 
  rows_per_strata <- c()
  for (k in 1:floodzip_id){
    rows_per_strata[k] = durations[k] * case_control_set #obtain number of rows in each strata 
    N = N + rows_per_strata[k]
  } 
  
  #obtain starting point for each new strata 
  floodzip_id_dat <- dat %>% select(-c(2:11)) %>% mutate(index = row_number()) 
  sequence <- which(floodzip_id_dat$index == 1)
  
  offset <- dat$population
  
  Sigma_d <- matrix(NA, nrow = num_coeff, ncol = num_coeff)
  Sigma_t <- matrix(NA, nrow = num_coeff, ncol = num_coeff)
  for (i in 1:num_coeff) {
    for (j in 1:num_coeff) {
      d_i <- ceiling((sqrt(1 + 8*i) - 1)/2)
      d_j <- ceiling((sqrt(1 + 8*j) - 1)/2)
      t_i <- i - d_i*(d_i-1)/2
      t_j <- j - d_j*(d_j-1)/2
      Sigma_d[i,j] <- abs(d_i - d_j)
      Sigma_t[i,j] <- abs(t_i - t_j)
    }
  }
  
  Sigma = exp((-1/phi)*Sigma_d + (-1/tau)*Sigma_t + log(little_sigma2)) #multi-dimensional product
  
  X <- readRDS(paste0(path_to_data, '/X_testdata.rds'))
  
  mu <- rep(0, num_coeff)
  
  return(list(floodzip_id = floodzip_id,
              case_control_set = case_control_set,
              durations = durations,
              D = D,
              num_coeff = num_coeff,
              N = N,
              rows_per_strata = rows_per_strata,
              sequence = sequence,
              X = X,
              offset = offset,
              Sigma_d = Sigma_d,
              Sigma_t = Sigma_t,
              mu = mu,
              Sigma = Sigma))
}


simulate_beta <- function(mu, Sigma){
  
  library(mvtnorm)
  set.seed(2024)
  true_beta <- as.vector(mvtnorm::rmvnorm(n = 1, mean = mu, sigma = Sigma))
  
  return(true_beta)
}


get_prob <- function(N, X, beta, offset, sequence, strata, rows_per_strata){
  
  log_numer <- vector()
  
  for (obs in 1:N){
    log_numer[obs] = X[obs,] %*% beta + log(offset[obs])
  }
  
  for (obs in 1:N){
    log_numer[obs] = X[obs,] %*% true_beta + log(offset[obs])
  }
  
  log_denom <- vector()
  prob <- vector()
  
  library(matrixStats)
  
  for (strata in 1:length(sequence)){
    log_denom[strata] = logSumExp(log_numer[(sequence[strata]): (sequence[strata] + rows_per_strata[strata] - 1)])
    prob[(sequence[strata]): (sequence[strata] + rows_per_strata[strata] - 1)] = exp(log_numer[(sequence[strata]): (sequence[strata] + rows_per_strata[strata] - 1)] - log_denom[strata])
  }
  
  return(prob)
}

get_case_counts <- function(sequence, strata, rows_per_strata, prob){
  Y <- vector()
  
  sample_size <- sample((1:50), length(sequence), replace = TRUE)
  
  for (strata in 1:length(sequence)){
    Y[(sequence[strata]): (sequence[strata] + rows_per_strata[strata] - 1)] <- rmultinom(n = 1, 
                                                                                         size = sample_size[strata], 
                                                                                         prob = prob[(sequence[strata]): (sequence[strata] + rows_per_strata[strata] - 1)])
  }
  return(Y)
}

simulation <- 7
data_dir <- paste0('/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/data/simulations/general_simulation', simulation, '/')
ifelse(!dir.exists(data_dir), dir.create(data_dir, recursive=TRUE), FALSE)

data <- data_setup(path_to_data = '/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/data', 
                   little_sigma2 = 0.2,
                   phi = 0.5, 
                   tau = 0.5)

true_beta <- simulate_beta(mu = data$mu,
                           Sigma = data$Sigma)

data$Sigma <- NULL 
saveRDS(data, paste0(data_dir, 'allfloodzips_dur3_mult_gp_general_simulation', simulation, '.rds'))

saveRDS(true_beta, paste0(data_dir, 'allfloodzips_dur3_true_beta_general_simulation', simulation, '.rds'))

prob <- get_prob(N = data$N,
                 X = data$X,
                 beta = true_beta, 
                 offset = data$offset,
                 sequence = data$sequence,
                 strata = data$strata,
                 rows_per_strata = data$rows_per_strata)

nsim <- 5000
Y_mat <- matrix(data = NA, nrow = data$N, ncol = nsim)
for (i in 1:nsim){
  Y_mat[,i] <- get_case_counts(sequence = data$sequence, 
                               strata = data$strata, 
                               rows_per_strata = data$rows_per_strata,
                               prob = prob)
}

saveRDS(Y_mat, paste0(data_dir, '/allfloodzips_dur3_Y_general_simulation', simulation,'.rds'))



