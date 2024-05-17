# April 2024
# Generate data from the true DGP (one realization of the GP prior)

setwd('/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/') #cluster 

suppressMessages(library(tidyverse))

data_setup <- function(path_to_data, little_sigma_beta2, phi, tau, little_sigma_theta2, eta){
  dat <- readRDS(paste0(path_to_data, '/lags_test_data.rds'))
  dat <- dat %>% arrange(floodzip_id,zipcode,month,day,control_indicator) #order to match paper 
  floodzip_id <- length(unique(dat$floodzip_id)) %>% as.double()
  case_control_set <- length(unique(dat$control_indicator)) %>% as.double()
  durations <- dat %>% group_by(floodzip_id) %>% summarise(durations = unique(d)) %>% select(durations) %>% unname() %>% unlist()
  D <- max(durations) 
  num_coeff <- D*(D+1)/2
  lag <- 5
  durations_with_lag <- durations + lag 
  
  
  #define N as 0, then compute since there are multiple sets of rows per flood-zip combination if duration > 1
  #N <- nrow(dat)
  N <- 0 
  rows_per_strata <- c()
  for (k in 1:floodzip_id){
    rows_per_strata[k] = durations_with_lag[k] * case_control_set #obtain number of rows in each strata which now include lag days 
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
  
  Sigma_beta = exp((-1/phi)*Sigma_d + (-1/tau)*Sigma_t + log(little_sigma_beta2)) #multi-dimensional product
  
  Sigma_l <- matrix(NA, nrow = lag, ncol = lag)
  for (i in 1:lag){
    for (j in 1:lag){
      Sigma_l[i,j] <- abs(i - j)
    }
  }
  
  Sigma_theta = exp((-1/eta)*Sigma_l + log(little_sigma_theta2))
  
  #create X matrix: define as 0, then figure out which columns exposure is in based on values of "d" and "t"
  #note: only every 3rd row is an event
  
  # X <- matrix(data = 0, nrow = nrow(dat), ncol = num_coeff)
  # exposure_mapped <- apply(dat, 1, function(i) {t <- as.numeric(i[10]); d <- as.numeric(i[11]); d*(d-1)/2+t})[seq(1,N, by = 3)]
  # sequence_events <- seq(1,N, by = 3)
  # for (i in 1:length(sequence_events)){
  #   X[sequence_events[i],exposure_mapped[i]] <- 1
  # }
  # saveRDS(X, "/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/data/X_lags_testdata.rds")
  # 
  # A <- as.matrix(dat[,c("lag1", "lag2", "lag3", "lag4", "lag5")])
  # saveRDS(A, "/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/data/A_lags_testdata.rds")
  
  X <- readRDS(paste0(path_to_data, '/X_lags_testdata.rds'))
  A <- readRDS(paste0(path_to_data, '/A_lags_testdata.rds'))
  
  mu_beta <- rep(0, num_coeff)
  mu_theta <- rep(0, lag)
  
  return(list(floodzip_id = floodzip_id,
              case_control_set = case_control_set,
              durations = durations,
              D = D,
              num_coeff = num_coeff,
              lag = lag,
              N = N,
              rows_per_strata = rows_per_strata,
              sequence = sequence,
              X = X,
              A = A,
              offset = offset,
              Sigma_d = Sigma_d,
              Sigma_t = Sigma_t,
              Sigma_l = Sigma_l,
              mu_beta = mu_beta,
              mu_theta = mu_theta,
              Sigma_beta = Sigma_beta,
              Sigma_theta = Sigma_theta))
}


simulate_beta <- function(mu, Sigma){
  
  library(mvtnorm)
  set.seed(2024)
  true_beta <- as.vector(mvtnorm::rmvnorm(n = 1, mean = mu, sigma = Sigma))
  
  return(true_beta)
}

simulate_theta <- function(mu, Sigma){
  library(mvtnorm)
  set.seed(2024)
  true_theta <- as.vector(mvtnorm::rmvnorm(n = 1, mean = mu, sigma = Sigma))
  
  return(true_theta)
}


get_prob <- function(N, X, A, beta, theta, offset, sequence, strata, rows_per_strata){
  
  log_numer <- vector()
  
  for (obs in 1:N){
    log_numer[obs] = X[obs,] %*% beta + A[obs,] %*% theta + log(offset[obs])
  }
  
  # don't think these lines get utilized 
  # for (obs in 1:N){
  #   log_numer[obs] = X[obs,] %*% true_beta + A[obs, ] %*% true_theta + log(offset[obs])
  # }
  
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

simulation <- 14
data_dir <- paste0('/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/data/simulations/general_simulation', simulation, '/')
ifelse(!dir.exists(data_dir), dir.create(data_dir, recursive=TRUE), FALSE)

data <- data_setup(path_to_data = '/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/data', 
                   little_sigma_beta2 = 0.2,
                   phi = 0.4, 
                   tau = 0.6,
                   little_sigma_theta2 = 0.1,
                   eta = 0.5)

true_beta <- simulate_beta(mu = data$mu_beta,
                           Sigma = data$Sigma_beta)

true_theta <- simulate_theta(mu = data$mu_theta,
                             Sigma = data$Sigma_theta)

data$Sigma_beta <- NULL 
data$Sigma_theta <- NULL 

saveRDS(data, paste0(data_dir, 'allfloodzips_dur3_lags_mult_gp_general_simulation', simulation, '.rds'))

saveRDS(true_beta, paste0(data_dir, 'allfloodzips_dur3_lags_true_beta_general_simulation', simulation, '.rds'))
saveRDS(true_theta, paste0(data_dir, 'allfloodzips_dur3_lags_true_theta_general_simulation', simulation, '.rds'))

prob <- get_prob(N = data$N,
                 X = data$X,
                 A = data$A, 
                 beta = true_beta, 
                 theta = true_theta,
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

saveRDS(Y_mat, paste0(data_dir, '/allfloodzips_dur3_lags_Y_general_simulation', simulation,'.rds'))


