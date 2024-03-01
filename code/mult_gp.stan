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
  int<lower=1> sequence[N / case_control_set]; // generate a sequence of numbers to indicate the start of a new strata (i.e. every three rows)
  
  vector[num_coeff] mu; // mu is a vector of length D*(D+1)/2
}

transformed data {
  matrix[num_coeff, num_coeff] log_Sigma; // set as fixed to understand where issue is
  log_Sigma = -Sigma_d - Sigma_t;
}

// The parameters accepted by the model.
parameters {
  real little_sigma; // non-zero 
  // real phi; // non-zero, lengthscale for d
  // real tau; // non-zero, lengthscale for t
  vector[num_coeff] beta; // beta is a vector of length D*(D+1)/2
}

transformed parameters {
  // note: 'numer', 'denom', and 'prob' get sampled (could be slower); 
  // vector[N] numer; 
  // vector[N / case_control_set] denom;
  
  vector[N] log_numer;
  vector[N / case_control_set] log_denom;
  simplex[case_control_set] prob[N/case_control_set]; 
  
  //for (obs in 1:N){
    //numer[obs] = exp(X[obs] * beta) * offset[obs]; 
  //}
  
  //for (strata in 1:num_elements(sequence)){
    //RACHEL :: not sure about this but, instead of '+2' in the rows below, for the general case shouldn't it be '+case_control_set'
    //denom[strata] = sum(numer[(sequence[strata]): (sequence[strata]+ 2)]); // compute denominator for each observation in a strata only once
    //prob[strata] = numer[(sequence[strata]): (sequence[strata] + 2)]/denom[strata]; // divide for each observation and save probabilities in a column of matrix
  //}
  
  for (obs in 1:N){
    log_numer[obs] = X[obs] * beta + log(offset[obs]); 
  }
  
  for (strata in 1:num_elements(sequence)){
    // RACHEL :: not sure about this but, instead of '+2' in the rows below, for the general case shouldn't it be '+case_control_set'
    log_denom[strata] = log_sum_exp(log_numer[(sequence[strata]): (sequence[strata] + 2)]); // compute denominator for each observation in a strata only once
    prob[strata] = exp(log_numer[(sequence[strata]): (sequence[strata] + 2)] - log_denom[strata]); // divide for each observation and save probabilities in a column of matrix
  }
  
  // matrix[num_coeff, num_coeff] Sigma = little_sigma*exp(-1/phi*Sigma_d + (-1/tau*Sigma_t)) where Sigma_d[i,j] = |d_i - d_j|, Sigma_t[i,j] = |t_i - t_j|
  matrix[num_coeff, num_coeff] log_Sigma2 = log_Sigma + log(little_sigma);
}
  
model {
  // likelihood 
  for (strata in 1:num_elements(sequence)){
    target += multinomial_lupmf(Y[(sequence[strata]): (sequence[strata] + 2)] | prob[strata]); // sum over log pmf by strata up to additive constants
  }

  // priors
  little_sigma ~ inv_gamma(5,5); // parameters for covariance cause Sigma to not be symmetric (?) 
  // phi ~ inv_gamma(5,5);
  // tau ~ inv_gamma(5,5);
  beta ~ multi_normal(mu, exp(log_Sigma2)); // beta is a vector of length D*(D+1)/2
}
