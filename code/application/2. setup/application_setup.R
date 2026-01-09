# April 2024

setwd('./multinomial_GP/') #cluster 

suppressMessages(library(tidyverse))

data_setup <- function(path_to_data){
  dat <- readRDS(paste0(path_to_data, '/medicare_cause_2000_2016_no_lag_dur10.rds')) #same in 1d and 2d
  dat <- dat %>% arrange(floodcty_id,fipscounty,month,day,control_indicator) #order to match paper 
  floodcty_id <- length(unique(dat$floodcty_id)) %>% as.double()
  case_control_set <- length(unique(dat$control_indicator)) %>% as.double()
  durations <- dat %>% group_by(floodcty_id) %>% summarise(durations = unique(duration)) %>% select(durations) %>% unname() %>% unlist()
  D <- max(durations) 
  num_coeff <- D*(D+1)/2
  #lag <- 5
  #num_lag_coeff <- lag #this implies that we estimate coefficients for all durations in spite of missing ones 
  num_conf <- 6 #number of confounders
  spline_df <- 3 #degrees of freedom for spline 
  #durations_with_lag <- durations + lag 
  
  
  #define N as 0, then compute since there are multiple sets of rows per flood-zip combination if duration > 1
  #N <- nrow(dat)
  N <- 0 
  rows_per_strata <- c()
  for (k in 1:floodcty_id){
    #rows_per_strata[k] = durations_with_lag[k] * case_control_set #obtain number of rows in each strata which now include lag days (change to durations_with_lag[k]) 
    rows_per_strata[k] = durations[k] * case_control_set #obtain number of rows in each strata which now include lag days (change to durations_with_lag[k]) 
    N = N + rows_per_strata[k]
  } 
  
  #obtain starting point for each new strata 
  sequence <- as.numeric(rownames(dat[match(unique(dat$floodcty_id), dat$floodcty_id),]))
  
  cases <- dat$cases
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
  # Sigma_d2 <- matrix(NA, nrow = num_lag_coeff, ncol = num_lag_coeff)
  # Sigma_l <- matrix(NA, nrow = num_lag_coeff, ncol = num_lag_coeff)
  # for (i in 1:num_lag_coeff){
  #   for (j in 1:num_lag_coeff){
  #     Sigma_l[i,j] <- abs(i - j)
  #     d_i <- ceiling(i/lag)
  #     d_j <- ceiling(j/lag)
  #     l_i <- i - lag*(d_i - 1)
  #     l_j <- j - lag*(d_j - 1)
  #     Sigma_d2[i,j] <- abs(d_i - d_j)
  #     Sigma_l[i,j] <- abs(l_i - l_j)
  #   }
  # }
  
  #phi replaces gamma if same hyperparameter for duration for flooded days and lag days 
  
  #create X matrix: define as 0, then figure out which columns exposure is in based on values of "d" and "t"
  #note: only every 3rd row is an event
  
  #same in 1d and 2d
  X <- matrix(data = 0, nrow = nrow(dat), ncol = num_coeff)
  exposure_mapped <- apply(dat, 1, function(i) {t <- as.numeric(i[11]); d <- as.numeric(i[12]); d*(d-1)/2+t})[seq(1,N, by = 3)] #index might have to change depending on number of columns in dataset
  sequence_events <- seq(1,N, by = 3)
  for (i in 1:length(sequence_events)){
    X[sequence_events[i],exposure_mapped[i]] <- 1
  }
  #saveRDS(X, "X_cause_2000_2016_cov.rds"))
  
  #1d lags 
  #A <- as.matrix(dat[,c("lag1", "lag2", "lag3", "lag4", "lag5", "lag6", "lag7")])
  #saveRDS(A, "A_cause_cov_lags_1d.rds")
  
  #Z <- as.matrix(dat[,c("temp", "hum", "wind", "pm25", "no2", "ozone")]) #zipcode variable names
  #Z <- as.matrix(dat[,c("tmmx", "rmax", "vs", "pm25", "o3", "no2")]) #county variable names
  Z <- readRDS('./multinomial_GP/FFS_final_oct2025/spline_covariates_dur10.rds')
  
  #2d lags 
  # A <- matrix(data = 0, nrow = nrow(dat), ncol = num_lag_coeff)
  # lag_mapped <- apply(dat, 1, function(i) {l <- as.numeric(i[12]); d <- as.numeric(i[11]); lag*(d-1) + l})[seq(1,N,by = 3)]
  # sequence_lags <- seq(1,N,by=3)
  # for (i in 1:length(sequence_events)){
  #   A[sequence_lags[i],lag_mapped[i]] <- 1
  # }
  # saveRDS(A, "A_cause_cov_lags_2d.rds")

  mu_beta <- rep(0, num_coeff)
  #mu_theta <- rep(0, lag) #1d
  #mu_theta <- rep(0, num_lag_coeff) #2d
  
  #No LAGS 
  return(list(floodcty_id = floodcty_id,
              case_control_set = case_control_set,
              durations = durations,
              D = D,
              num_coeff = num_coeff,
              num_conf = num_conf,
              spline_df = spline_df,
              N = N,
              rows_per_strata = rows_per_strata,
              sequence = sequence,
              X = X,
              Z = Z,
              Y = cases,
              offset = offset,
              Sigma_d = Sigma_d,
              Sigma_t = Sigma_t,
              mu_beta = mu_beta))
  
  #1d LAGS 
  # return(list(floodcty_id = floodcty_id,
  #             case_control_set = case_control_set,
  #             durations = durations,
  #             D = D,
  #             num_coeff = num_coeff,
  #             lag = lag,
  #             num_lag_coeff = num_lag_coeff,
  #             num_conf = num_conf,
  #             spline_df = spline_df,
  #             N = N,
  #             rows_per_strata = rows_per_strata,
  #             sequence = sequence,
  #             X = X,
  #             A = A,
  #             Z = Z,
  #             Y = cases,
  #             offset = offset,
  #             Sigma_d = Sigma_d,
  #             Sigma_t = Sigma_t,
  #             Sigma_l = Sigma_l,
  #             mu_beta = mu_beta,
  #             mu_theta = mu_theta))
  
  #2d LAGS
  #   return(list(floodcty_id = floodcty_id,
  #               case_control_set = case_control_set,
  #               durations = durations,
  #               D = D,
  #               num_coeff = num_coeff,
  #               lag = lag,
  #               num_lag_coeff = num_lag_coeff,
  #               N = N,
  #               rows_per_strata = rows_per_strata,
  #               sequence = sequence,
  #               X = X,
  #               A = A,
  #               offset = offset,
  #               Sigma_d = Sigma_d,
  #               Sigma_t = Sigma_t,
  #               Sigma_d2 = Sigma_d2,
  #               Sigma_l = Sigma_l,
  #               mu_beta = mu_beta,
  #               mu_theta = mu_theta,
  #               Sigma_beta = Sigma_beta,
  #               Sigma_theta = Sigma_theta))
}

data_dir <- './multinomial_GP/FFS_final_oct2025/'
#ifelse(!dir.exists(data_dir), dir.create(data_dir, recursive=TRUE), FALSE)

data <- data_setup(path_to_data = data_dir)

saveRDS(data, paste0(data_dir, 'cause_2000_2016_spline_covariates_no_lag_dur10_mult_gp_real_data.rds'))

