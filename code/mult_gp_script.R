setwd('/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/') #cluster 

suppressMessages(library(tidyverse))

dat <- readRDS('data/test_data.rds')
dat <- dat %>% arrange(floodzip_id,zipcode,month,day,control_indicator) #order to match paper 

library(rstan)

floodzip_id <- length(unique(dat$floodzip_id)) %>% as.double()
case_control_set <- length(unique(dat$control_indicator)) %>% as.double()

durations <- dat %>% group_by(floodzip_id) %>% summarise(durations = unique(d)) %>% select(durations) %>% unname() %>% unlist()
D <- max(durations)
num_coeff <- D*(D+1)/2

#define N as 0, then compute since there are multiple sets of rows per flood-zip combination if duration > 1
N <- 0 
for (k in 1:floodzip_id){
  count = durations[k] * case_control_set
  N = N + count
} 

Y <- dat$cases

#create X matrix: define as 0, then figure out which columns exposure is in based on values of "d" and "t"
#note: only every 3rd row is an event

# X <- matrix(data = 0, nrow = nrow(dat), ncol = num_coeff)
# exposure_mapped <- apply(dat, 1, function(i) {t <- as.numeric(i[10]); d <- as.numeric(i[11]); d*(d-1)/2+t})[seq(1,N, by = 3)]
# sequence <- seq(1,N, by = 3)
# for (i in 1:length(sequence)){
#   X[sequence[i],exposure_mapped[i]] <- 1
# }
# saveRDS(X, "/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/data/X_testdata_175floodzips.rds")

X <- readRDS('data/X_testdata.rds')

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

sequence <- seq(1, N, by = 3)

mu <- rep(0, num_coeff)

mult_gp_data <- list(floodzip_id = floodzip_id,
                     case_control_set = case_control_set,
                     durations = durations,
                     D = D,
                     num_coeff = num_coeff,
                     N = N,
                     Y = Y,
                     X = X,
                     offset = offset,
                     Sigma_d = Sigma_d,
                     Sigma_t = Sigma_t,
                     sequence = sequence,
                     mu = mu)


options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

suppressMessages(
test <- stan(
  file = 'code/mult_gp.stan',  # Stan program
  data = mult_gp_data,    # named list of data
  chains = 4, 
  warmup = 500,
  iter = 2000
))

saveRDS(test, '/n/dominici_nsaph_l3/Lab/projects/floods-hospitalizations-glm/multinomial_GP/output/run.rds')
  
print(test, pars = c("beta", "little_sigma2", "phi", "tau", "lp__"))
pairs(test, pars = c("beta", "little_sigma2", "phi", "tau", "lp__"))







