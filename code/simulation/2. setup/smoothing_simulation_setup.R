# August 2024
# Generate data from a smoothing spline (number of basis functions controls smoothing)

setwd('/n/holylabs/nethery_lab/Lab/floods-hospitalizations-glm/multinomial_GP/') #cluster 

suppressMessages(library(tidyverse))

data_setup <- function(path_to_data, little_sigma2 = NULL, phi = NULL, tau = NULL){
  dat <- readRDS(paste0(path_to_data, '/rand_data.rds'))
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
  
  #create X matrix: define as 0, then figure out which columns exposure is in based on values of "d" and "t"
  #note: only every 3rd row is an event
  
  # X <- matrix(data = 0, nrow = nrow(dat), ncol = num_coeff)
  # exposure_mapped <- apply(dat, 1, function(i) {t <- as.numeric(i[10]); d <- as.numeric(i[11]); d*(d-1)/2+t})[seq(1,N, by = 3)]
  # sequence_events <- seq(1,N, by = 3)
  # for (i in 1:length(sequence_events)){
  #   X[sequence_events[i],exposure_mapped[i]] <- 1
  # }
  # saveRDS(X, 'X_randdata.rds')
  
  X <- readRDS(paste0(path_to_data, '/X_randdata.rds'))
  
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
              mu = mu))
}


smooth_beta <- function(num_bf, D, num_coeff){
   ntimes <- 1:D
   d <- rep(1:D, ntimes)
   t <- c(1, 1:(D-12), 1:(D-11), 1:(D-10), 1:(D-9), 1:(D-8), 1:(D-7), 1:(D-6), 1:(D-5), 1:(D-4), 1:(D-3), 1:(D-2), 1:(D-1), 1:D)

   coord <- matrix(c(d, t), byrow = F, ncol = 2)
   colnames(coord) <- c("d", "t")

   ##https://asbates.rbind.io/2019/03/01/thin-plate-splines/
   beta <- rep(1,num_coeff) #does not matter, just need something for the dummy model
   beta <- rnorm(n=num_coeff, mean=0, sd =1)
   coord_beta <- data.frame(coord, beta)

   library(mgcv)

   fit <- gam(coord_beta$beta ~ s(coord_beta$d, coord_beta$t, bs = "tp", k = num_bf + 1), fit=FALSE)

   X_2d <- fit$X[,-1] #remove intercept

   #generate magnitude for each basis function from std. normal
   coef <- rnorm(num_bf) #number of basis functions should equal the number of columns in the model matrix (+1 with intercept)

   #smooth function is the linear combination of the natural spline values and the magnitudes/coefficients
   smooth_fn <- X_2d %*% coef
   plot(smooth_fn)

   return(smooth_fn)
}


get_prob <- function(N, X, beta, offset, sequence, strata, rows_per_strata){
  
  log_numer <- vector()
  
  for (obs in 1:N){
    log_numer[obs] = X[obs,] %*% beta + log(offset[obs])
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
data_dir <- paste0('/n/holylabs/LABS/nethery_lab/Lab/floods-hospitalizations-glm/multinomial_GP/data/simulations/general_simulation', simulation, '/')
ifelse(!dir.exists(data_dir), dir.create(data_dir, recursive=TRUE), FALSE)

data <- data_setup(path_to_data = '/n/holylabs/LABS/nethery_lab/Lab/floods-hospitalizations-glm/multinomial_GP/data')

true_beta <- smooth_beta(num_bf = 5,
                         D = data$D,
                         num_coeff = data$num_coeff)

saveRDS(data, paste0(data_dir, 'rand_mult_gp_general_smooth_med_', simulation, '.rds'))

saveRDS(true_beta, paste0(data_dir, 'rand_true_beta_general_smooth_med_', simulation, '.rds'))

#add noise example
#true_beta_Xpct <- true_beta + rnorm(105, 0, sqrt((X/100)*var(true_beta))) #X is the percent noise added



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

saveRDS(Y_mat, paste0(data_dir, '/rand_Y_general_smooth_med_', simulation,'.rds'))



