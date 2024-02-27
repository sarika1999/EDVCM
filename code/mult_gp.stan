// February 2023
// Multinomial (conditional poisson) with GP prior with mean 'mu' and
// standard deviation 'Sigma' (sum of exponential kernels in two dimension).

data {
  int<lower=1> floodzip_id; // number of unique flood-zip combinations 
  int<lower=1> case_control_set; // number of cases (2) + controls (1) in a matched set
  vector<lower=1>[floodzip_id] durations; // duration for each unique flood-zip combination 
  int<lower=1> D; // max duration 
  int<lower=1> num_coeff; // define because it is used throughout the program 
  int<lower=1> N; // total number of rows in the dataset

  int<lower=0> Y[N]; // Y (cases) is a vector of length N 
  matrix<lower=0, upper = 1>[N,num_coeff] X; // X is a matrix with dimensions N x D*(D+1)/2, takes on values 0 or 1
  int<lower=1> offset[N]; // offset (population) is a vector of length N
  
  matrix[num_coeff,num_coeff] Sigma_d;
  matrix[num_coeff,num_coeff] Sigma_t;
  
  // N should be exactly divisible by case_control_set to produce an integer
  // int<lower=0> sum_Y[N / case_control_set]; // obtain sum of cases for every strata (i.e three rows)
  int<lower=1> sequence[N / case_control_set]; // generate a sequence of numbers to indicate the start of a new strata (i.e. every three rows)
  
  vector[num_coeff] mu; // mu is a vector of length D*(D+1)/2
}

// The parameters accepted by the model.
parameters {
  // real eta; // non-zero 
  // real phi; // non-zero, lengthscale for d
  // real tau; // non-zero, lengthscale for t
  vector[num_coeff] beta; // beta is a vector of length D*(D+1)/2
}

//transformed parameters {
  // real little_sigma = 1/eta; 
// }
  
model {
  // likelihood 
  vector[N] numer; 
  vector[N / case_control_set] denom;
  matrix[N / case_control_set, case_control_set] prob;
  // vector[N] lik;
  matrix[num_coeff, num_coeff] log_Sigma = -Sigma_d - Sigma_t; // set as fixed to understand where issue is
  // matrix[num_coeff, num_coeff] Sigma = little_sigma*exp(-1/phi*Sigma_d + (-1/tau*Sigma_t)) where Sigma_d[i,j] = |d_i - d_j|, Sigma_t[i,j] = |t_i - t_j|
  // matrix[num_coeff, num_coeff] log_Sigma = log(little_sigma) - 1/phi*Sigma_d - 1/tau*Sigma_t;
  
  for (obs in 1:N){
    numer[obs] = exp(X[obs] * beta) * offset[obs]; 
  }
  for (strata in 1:num_elements(sequence)){
    denom[strata] = sum(numer[(sequence[strata]): (sequence[strata] + 2)]); // compute denominator for each observation in a strata only once 
    prob[strata,] = to_row_vector(numer[(sequence[strata]): (sequence[strata] + 2)]/denom[strata]); // divide for each observation and save probabilities in a row of matrix
    // Y[(sequence[strata]): (sequence[strata] + 2)] ~ multinomial(to_vector(prob[strata,])); // needs to be a product for ALL data 
    // multinomial_lupmf() is not available in 2.21
    target += multinomial_lpmf(Y[(sequence[strata]): (sequence[strata] + 2)] | to_vector(prob[strata,])); // sum over log pmf for each strata 
  }
  
  // priors
  // eta ~ inv_gamma(5, 5); 
  // phi ~ inv_gamma(5, 5);
  // tau ~ inv)gamma(5, 5);
  beta ~ multi_normal(mu, exp(log_Sigma)); // beta is a vector of length D*(D+1)/2
}

