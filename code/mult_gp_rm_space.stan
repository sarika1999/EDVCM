data {
  int<lower=1> floodzip_id; // number of unique flood-zip combinations 
  int<lower=1> case_control_set; // number of cases (2) + controls (1) in a matched set
  vector<lower=1>[floodzip_id] durations; // duration for each unique flood-zip combination 
  int<lower=1> D; // max duration 
  int<lower=1> num_coeff; // define because it is used throughout the program 
  int<lower=1> N; // compute in R script and pass into STAN 
  int<lower=0> Y[N]; // Y (cases) is a vector of length N 
  matrix<lower=0, upper = 1>[N,num_coeff] X; // X is a matrix with dimensions N x D*(D+1)/2
  int<lower=1> offset[N]; // offset (population) is a vector of length N
  matrix[num_coeff*num_coeff, 4] td_combos; // possible combinations of time-point and duration values within a flood event 
  int<lower=0> sum_Y[N / case_control_set]; // obtain sum of cases for every strata (i.e three rows)
  int<lower=1> sequence[N / case_control_set]; // generate a sequence of numbers to indicate the start of a new strata (i.e. every three rows)
  vector<lower=0>[num_coeff] mu; // mu is a vector of length D*(D+1)/2
  real sigma_t; // sigma for t == 1
  real sigma_d; // sigma for d == 1
}
transformed data {
  real d_i[num_coeff*num_coeff] = to_array_1d(td_combos[,1]);  
  real t_i[num_coeff*num_coeff] = to_array_1d(td_combos[,2]);
  real d_j[num_coeff*num_coeff] = to_array_1d(td_combos[,3]);
  real t_j[num_coeff*num_coeff] = to_array_1d(td_combos[,4]);
}
parameters {
  real<lower=1e-9> eta; //non-zero
  real<lower=1e-9> phi; // non-zero, lengthscale for d
  real<lower=1e-9> tau; // non-zero, lengthscale for t
  vector[num_coeff] beta; // beta is a vector of length D*(D+1)/2
}
transformed parameters {
  real little_sigma = 1/eta; 
}
model {
  // likelihood 
  vector[N] numer; 
  vector[N] denom;
  vector[N] prob;
  vector[N] log_lik;
  matrix[num_coeff, num_coeff] Sigma = little_sigma*gp_exponential_cov(d_i, d_j, sigma_d, phi)*gp_exponential_cov(t_i, t_j, sigma_t, tau); 
  for (strata in 1:num_elements(sequence)){ //no length in STAN 
    for (obs in (sequence[strata]): (sequence[strata] + 2)){
      // j = 1,...,N indexes flood-zipcode-day
      // Y ~ Poisson(lambda_j) --> log(lambda_j) = X_j * beta + log(offset_j)
      numer[obs] = exp(X[obs] * beta) * offset[obs]; 
      // log(pi_j) = X_j * beta + log(offset_j) - log(sum_k(exp(X_k * beta)*offset_k))
      denom[obs] = sum(exp(numer[(sequence[strata]): (sequence[strata] + 2)]));
      // pi_j = exp(X_j * beta)*offset_j / sum_k(exp(X_k * beta)*offset_k)
      prob[obs] = numer[obs]/denom[obs];
      log_lik[obs] = tgamma(sum_Y[strata])^(1 / case_control_set)*(prob[obs])^(Y[obs])*inv(tgamma(Y[obs]));
    }
  }
  eta ~ gamma(0.001, 0.001);
  phi ~ gamma(0.001, 0.001);
  tau ~ gamma(0.001, 0.001);
  beta ~ multi_normal(mu, Sigma); // beta is a vector of length D*(D+1)/2
}

