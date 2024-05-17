// April,May 2024
// Multinomial (conditional poisson) with GP prior with mean 'mu_beta' and
// standard deviation 'Sigma_beta' (product of exponential kernels in two dimension).
// with post-flood time (lags) that has a GP prior with mean 'mu_theta' and 
// standard deviation 'Sigma_theta' (exponential kernel in one dimension)

data {
  int<lower=1> floodzip_id; // number of unique flood-zip combinations 
  int<lower=1> case_control_set; // number of cases (2) + controls (1) in a matched set
  vector<lower=1>[floodzip_id] durations; // duration for each unique flood-zip combination 
  int<lower=1> D; // max duration 
  int<lower=1> num_coeff; // define because it is used throughout the program
  int<lower=1> lag; // number of lag days -- assume they are not a function of duration 
  int<lower=1> N; // total number of rows in the dataset
  
  int<lower=1> rows_per_strata[floodzip_id]; // case_control_set*durations[floodzip_id]
  int<lower=0> Y[N]; // Y (cases) is a vector of length N 
  matrix<lower=0, upper = 1>[N,num_coeff] X; // X is a matrix with dimensions N x D*(D+1)/2, takes on values 0 or 1
  matrix<lower=0, upper = 1>[N,lag] A; // A is a matrix with dimensions N x lag, takes on values 0 or 1
  int<lower=1> offset[N]; // offset (population) is a vector of length N
  
  
  matrix[num_coeff,num_coeff] Sigma_d;
  matrix[num_coeff,num_coeff] Sigma_t;
  matrix[lag,lag] Sigma_l;
  
  // N should be exactly divisible by case_control_set to produce an integer
  int<lower=1> sequence[floodzip_id]; // generate a sequence of numbers to indicate the start of a new strata
  
  vector[num_coeff] mu_beta; // mu_beta is a vector of length D*(D+1)/2
  vector[lag] mu_theta; // mu_theta is a vector of length lag 
}

// The parameters accepted by the model.
parameters {
  real little_sigma_beta2; // non-zero 
  real phi; // non-zero, lengthscale for d
  real tau; // non-zero, lengthscale for t
  
  real little_sigma_theta2; // non-zero
  real eta; // non-zero, lengthscale for l 
  
  vector[num_coeff] beta; // beta is a vector of length D*(D+1)/2
  vector[lag] theta; // theta is a vector of length lag (arbitrary number of days added)
}
  
model {
  matrix[num_coeff, num_coeff] Sigma_beta; // cannot define cov_matrix in model block
  matrix[lag, lag] Sigma_theta; 
  vector[N] log_numer;
  vector[floodzip_id] log_denom;
  vector[N] prob;
  
  // matrix[num_coeff, num_coeff] Sigma = little_sigma*exp(-1/phi*Sigma_d + (-1/tau*Sigma_t)); where Sigma_d[i,j] = |d_i - d_j|, Sigma_t[i,j] = |t_i - t_j|
  Sigma_beta = exp((-1/phi)*Sigma_d + (-1/tau)*Sigma_t + log(little_sigma_beta2)); // multi-dimensional product 
  Sigma_theta = exp((-1/eta)*Sigma_l + log(little_sigma_theta2)); // one-dimensional 
  
  for (obs in 1:N){
    log_numer[obs] = X[obs] * beta + A[obs] * theta + log(offset[obs]); 
  }
  
  for (strata in 1:(num_elements(sequence))){
    log_denom[strata] = log_sum_exp(log_numer[(sequence[strata]): (sequence[strata] + rows_per_strata[strata] - 1)]); // compute denominator for each observation in a strata only once
    prob[(sequence[strata]): (sequence[strata] + rows_per_strata[strata] - 1)] = exp(log_numer[(sequence[strata]): (sequence[strata] + rows_per_strata[strata] - 1)] - log_denom[strata]); // divide for each observation and save probabilities 
  }
  
  // likelihood 
  for (strata in 1:num_elements(sequence)){
    target += multinomial_lupmf(Y[(sequence[strata]): (sequence[strata] + rows_per_strata[strata] - 1)] | prob[(sequence[strata]): (sequence[strata] + rows_per_strata[strata] - 1)]); // sum over log pmf by strata up to additive constants
  }

  // priors
  little_sigma_beta2 ~ inv_gamma(5,5); // uniform(0,1); // lognormal(log(0.3), 0.2);
  phi ~ inv_gamma(5,5); // uniform(0,1); // lognormal(log(0.3), 0.2);
  tau ~ inv_gamma(5,5); // uniform(0,1) // lognormal(log(0.3), 0.2);
  little_sigma_theta2 ~ inv_gamma(5,5); // uniform(0,1); // lognormal(log(0.3), 0.2);
  eta ~ inv_gamma(5,5); // uniform(0,1) // lognormal(log(0.3), 0.2);
  beta ~ multi_normal(mu_beta, Sigma_beta); // beta is a vector of length D*(D+1)/2
  theta ~ multi_normal(mu_theta, Sigma_theta); // theta is a vector of length lag 
}
