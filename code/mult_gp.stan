// February 2023
// Multinomial (conditional poisson) with GP kernel
// with mean 'mu' and standard deviation 'Sigma'.

// 
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
  
  positive_ordered[D] t; // possible time-points within a flood event
  positive_ordered[D] d; // possible durations of a flood event 
  
  // N should be exactly divisible by case_control_set to produce an integer, but use integer division operator to appease STAN
  int<lower=0> sum_Y[N / case_control_set]; // obtain sum of cases for every strata (i.e three rows)
  int<lower=1> sequence[N / case_control_set]; // generate a sequence of numbers to indicate the start of a new strata (i.e. every three rows)
  
  vector<lower=0>[num_coeff] mu; // mu is a vector of length D*(D+1)/2
  
  real sigma_t; // sigma for t == 1
  real sigma_d; // sigma for d == 1
}

transformed data {
  real t_arr[D] = to_array_1d(t);  
  real d_arr[D] = to_array_1d(d); 
}

// The parameters accepted by the model.
parameters {
  real<lower=1e-9> eta; //non-zero
  real<lower=1e-9> phi; // non-zero, lengthscale for t
  real<lower=1e-9> tau; // non-zero, lengthscale for d
  vector[num_coeff] beta; // beta is a vector of length D*(D+1)/2
}

transformed parameters {
  real little_sigma = 1/eta; 
}
  
// The model to be estimated. 
model {
  // likelihood 
  vector[N] numer; 
  vector[N] denom;
  vector[N] prob;
  vector[N] log_lik;
  // if there are not consecutive durations 'd' (2, 3, 5) such that the dimensions match t' (1, 2, 3, 4, 5), then the covariance calculation will fail  
  matrix<lower=0>[num_coeff, num_coeff] Sigma = little_sigma*gp_exponential_cov(t_arr,sigma_t, phi)*gp_exponential_cov(d_arr,sigma_d, tau); 
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
  
  // priors
  eta ~ gamma(0.001, 0.001);
  // could change based on how correlated we want the daily effects to be (for example drop below 0.05 if 10-14 days apart)
  // idea from H. Chang, but unsure how to execute similarly
  phi ~ gamma(0.001, 0.001);
  tau ~ gamma(0.001, 0.001);
  beta ~ multi_normal(mu, Sigma); // beta is a vector of length D*(D+1)/2
}

