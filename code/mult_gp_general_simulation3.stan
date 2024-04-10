// April 2024
// Multinomial (conditional poisson) with GP prior with mean 'mu' and
// standard deviation 'Sigma' (product of exponential kernels in two dimension).

data {
  int<lower=1> floodzip_id; // number of unique flood-zip combinations 
  int<lower=1> case_control_set; // number of cases (2) + controls (1) in a matched set
  vector<lower=1>[floodzip_id] durations; // duration for each unique flood-zip combination 
  int<lower=1> D; // max duration 
  int<lower=1> num_coeff; // define because it is used throughout the program
  int<lower=1> N; // total number of rows in the dataset
  
  int<lower=1> rows_per_strata[floodzip_id]; // case_control_set*durations[floodzip_id]
  int<lower=0> Y[N]; // Y (cases) is a vector of length N 
  matrix<lower=0, upper = 1>[N,num_coeff] X; // X is a matrix with dimensions N x D*(D+1)/2, takes on values 0 or 1
  int<lower=1> offset[N]; // offset (population) is a vector of length N
  
  matrix[num_coeff,num_coeff] Sigma;
  
  // N should be exactly divisible by case_control_set to produce an integer
  int<lower=1> sequence[floodzip_id]; // generate a sequence of numbers to indicate the start of a new strata
  
  vector[num_coeff] mu; // mu is a vector of length D*(D+1)/2
}

// The parameters accepted by the model.
parameters {
  vector[num_coeff] beta; // beta is a vector of length D*(D+1)/2
}
  
model {
  vector[N] log_numer;
  vector[floodzip_id] log_denom;
  vector[N] prob;

  for (obs in 1:N){
    log_numer[obs] = X[obs] * beta + log(offset[obs]); 
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
  beta ~ multi_normal(mu, Sigma); // beta is a vector of length D*(D+1)/2
}
