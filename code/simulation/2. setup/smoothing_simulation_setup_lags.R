# April 2024
# Generate data from the true DGP (one realization of the GP prior)

setwd('/n/holylabs/nethery_lab/Lab/floods-hospitalizations-glm/multinomial_GP/') #cluster 

suppressMessages(library(tidyverse))

data_setup <- function(path_to_data, little_sigma_beta2 = NULL, phi = NULL, tau = NULL, little_sigma_theta2 = NULL, gamma = NULL, eta = NULL){
  dat <- readRDS(paste0(path_to_data, '/lags_rand_data.rds')) #same in 1d and 2d
  dat <- dat %>% arrange(floodzip_id,zipcode,month,day,control_indicator) #order to match paper 
  floodzip_id <- length(unique(dat$floodzip_id)) %>% as.double()
  case_control_set <- length(unique(dat$control_indicator)) %>% as.double()
  durations <- dat %>% group_by(floodzip_id) %>% summarise(durations = unique(d)) %>% select(durations) %>% unname() %>% unlist()
  D <- max(durations) 
  num_coeff <- D*(D+1)/2
  lag <- 5
  num_lag_coeff <- D*lag #this implies that we estimate coefficients for all durations in spite of missing ones 
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
  
  # 1d LAGS
  # Sigma_l <- matrix(NA, nrow = lag, ncol = lag)
  # for (i in 1:lag){
  #   for (j in 1:lag){
  #     Sigma_l[i,j] <- abs(i - j)
  #   }
  # }
  
  # 2d LAGS
  # Note: This is a different mapping function than main effects because we do not have a triangular coefficient structure for lags 
  Sigma_d2 <- matrix(NA, nrow = num_lag_coeff, ncol = num_lag_coeff)
  Sigma_l <- matrix(NA, nrow = num_lag_coeff, ncol = num_lag_coeff)
  for (i in 1:num_lag_coeff){
    for (j in 1:num_lag_coeff){
      Sigma_l[i,j] <- abs(i - j)
      d_i <- ceiling(i/lag)
      d_j <- ceiling(j/lag)
      l_i <- i - lag*(d_i - 1)
      l_j <- j - lag*(d_j - 1)
      Sigma_d2[i,j] <- abs(d_i - d_j)
      Sigma_l[i,j] <- abs(l_i - l_j)
    }
  }
  
  #phi replaces gamma if same hyperparameter for duration for flooded days and lag days 
  
  #create X matrix: define as 0, then figure out which columns exposure is in based on values of "d" and "t"
  #note: only every 3rd row is an event
  
  #same in 1d and 2d
  # X <- matrix(data = 0, nrow = nrow(dat), ncol = num_coeff)
  # exposure_mapped <- apply(dat, 1, function(i) {t <- as.numeric(i[10]); d <- as.numeric(i[11]); d*(d-1)/2+t})[seq(1,N, by = 3)]
  # sequence_events <- seq(1,N, by = 3)
  # for (i in 1:length(sequence_events)){
  #   X[sequence_events[i],exposure_mapped[i]] <- 1
  # }
  # saveRDS(X, "X_lags_2d_randdata.rds")
  
  #1d lags 
  # A <- as.matrix(dat[,c("lag1", "lag2", "lag3", "lag4", "lag5")])
  # saveRDS(A, "/A_lags_1d_randdata.rds")
  
  #2d lags 
  A <- matrix(data = 0, nrow = nrow(dat), ncol = num_lag_coeff)
  lag_mapped <- apply(dat, 1, function(i) {l <- as.numeric(i[12]); d <- as.numeric(i[11]); lag*(d-1) + l})[seq(1,N,by = 3)]
  sequence_lags <- seq(1,N,by=3)
  for (i in 1:length(sequence_events)){
    A[sequence_lags[i],lag_mapped[i]] <- 1
  }
  saveRDS(A, "A_lags_2d_randdata.rds")
  
  #X <- readRDS(paste0(path_to_data, '/X_lags_2d_randdata.rds'))
  
  mu_beta <- rep(0, num_coeff)
  #mu_theta <- rep(0, lag) #1d
  mu_theta <- rep(0, num_lag_coeff) #2d
  
  #1d LAGS 
  # return(list(floodzip_id = floodzip_id,
  #             case_control_set = case_control_set,
  #             durations = durations,
  #             D = D,
  #             num_coeff = num_coeff,
  #             lag = lag,
  #             N = N,
  #             rows_per_strata = rows_per_strata,
  #             sequence = sequence,
  #             X = X,
  #             A = A,
  #             offset = offset,
  #             Sigma_d = Sigma_d,
  #             Sigma_t = Sigma_t,
  #             Sigma_l = Sigma_l,
  #             mu_beta = mu_beta,
  #             mu_theta = mu_theta))
  
  #2d LAGS
  return(list(floodzip_id = floodzip_id,
              case_control_set = case_control_set,
              durations = durations,
              D = D,
              num_coeff = num_coeff,
              lag = lag,
              num_lag_coeff = num_lag_coeff,
              N = N,
              rows_per_strata = rows_per_strata,
              sequence = sequence,
              X = X,
              A = A,
              offset = offset,
              Sigma_d = Sigma_d,
              Sigma_t = Sigma_t,
              Sigma_d2 = Sigma_d2,
              Sigma_l = Sigma_l,
              mu_beta = mu_beta,
              mu_theta = mu_theta))
}


# use same betas as in smooth + 25% noise, smooth + 100% noise
# smooth_beta <- function(num_bf, D, num_coeff){
#   ntimes <- 1:D
#   d <- rep(1:D, ntimes)
#   t <- c(1, 1:(D-12), 1:(D-11), 1:(D-10), 1:(D-9), 1:(D-8), 1:(D-7), 1:(D-6), 1:(D-5), 1:(D-4), 1:(D-3), 1:(D-2), 1:(D-1), 1:D)
#   
#   coord <- matrix(c(d, t), byrow = F, ncol = 2)
#   colnames(coord) <- c("d", "t")
#   #
#   #   ##https://asbates.rbind.io/2019/03/01/thin-plate-splines/
#   beta <- rep(1,num_coeff) #does not matter, just need something for the dummy model
#   beta <- rnorm(n=num_coeff, mean=0, sd =1)
#   coord_beta <- data.frame(coord, beta)
#   #
#   library(mgcv)
#   #
#   fit <- gam(coord_beta$beta ~ s(coord_beta$d, coord_beta$t, bs = "tp", k = num_bf + 1), fit=FALSE)
#   #
#   X_2d <- fit$X[,-1] #remove intercept
#   #
#   #   #generate magnitude for each basis function from std. normal
#   coef <- rnorm(num_bf) #number of basis functions should equal the number of columns in the model matrix (+1 with intercept)
#   #
#   #   #smooth function is the linear combination of the natural spline values and the magnitudes/coefficients
#   smooth_fn <- X_2d %*% coef
#   plot(smooth_fn)
#   return(smooth_fn)
# }


smooth_theta <- function(num_bf, D, lag, num_lag_coeff){
  d <- rep(c(1:D), each =lag)
  t <- rep(1:lag, D)
  
  coord <- matrix(c(d, t), byrow = F, ncol = 2)
  colnames(coord) <- c("d", "t")
  #
  #   ##https://asbates.rbind.io/2019/03/01/thin-plate-splines/
  theta <- rep(1,num_lag_coeff) #does not matter, just need something for the dummy model
  theta <- rnorm(n=num_lag_coeff, mean=0, sd =1)
  coord_theta <- data.frame(coord, theta)
  #
  library(mgcv)
  #
  fit <- gam(coord_theta$theta ~ s(coord_theta$d, coord_theta$t, bs = "tp", k = num_bf + 1), fit=FALSE)
  #
  X_2d <- fit$X[,-1] #remove intercept
  #
  #   #generate magnitude for each basis function from std. normal
  coef <- rnorm(num_bf) #number of basis functions should equal the number of columns in the model matrix (+1 with intercept)
  #
  #   #smooth function is the linear combination of the natural spline values and the magnitudes/coefficients
  smooth_fn <- X_2d %*% coef
  plot(smooth_fn)
  
  smooth_fn_theta <- smooth_fn 
  return(smooth_fn_theta)
}

get_prob <- function(N, X, A, beta, theta, offset, sequence, strata, rows_per_strata){
  
  log_numer <- vector()
  
  for (obs in 1:N){
    log_numer[obs] = X[obs,] %*% beta + A[obs,] %*% theta + log(offset[obs])
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

simulation <- 
data_dir <- paste0('/n/holylabs/nethery_lab/Lab/floods-hospitalizations-glm/multinomial_GP/data/simulations/general_simulation', simulation, '/')
ifelse(!dir.exists(data_dir), dir.create(data_dir, recursive=TRUE), FALSE)

data <- data_setup(path_to_data = '/n/holylabs/nethery_lab/Lab/floods-hospitalizations-glm/multinomial_GP/data')

true_beta <- readRDS('./data/simulations/general_simulation/rand_true_beta_general_smooth_med_Xpctnoise_simulation.rds')
true_theta <- smooth_theta(num_bf = 5, D = 14, lag = 5, num_lag_coeff = D*lag)
#true_theta <- readRDS('./data/true_theta.rds')
true_theta <- true_theta + rnorm(D*lag, 0, sqrt((X/100)*var(true_theta))) #X is the percent noise added

saveRDS(data, paste0(data_dir, 'rand_mult_gp_general_smooth_med_simulation_Xpctnoise_2d_lags', simulation, '.rds'))

saveRDS(true_theta, paste0(data_dir, 'rand_true_theta_general_smooth_med_simulation_Xpctnoise_2d_lags', simulation, '.rds'))

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

saveRDS(Y_mat, paste0(data_dir, '/rand_Y_general_smooth_med_simulation_Xpctnoise_2d_lags', simulation,'.rds'))
