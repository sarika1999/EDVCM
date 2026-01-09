// Multinomial (conditional poisson) with GP prior with mean 'mu_beta' and
// standard deviation 'Sigma_beta' (product of exponential kernels in two dimension).
// with post-flood time (lags) that has a GP prior with mean 'mu_theta' and 
// standard deviation 'Sigma_theta' (product of exponential kernels in two dimension)

data {
  int<lower=1> floodzip_id; // number of unique flood-zip combinations 
  int<lower=1> case_control_set; // number of cases (2) + controls (1) in a matched set
  vector<lower=1>[floodzip_id] durations; // duration for each unique flood-zip combination 
  int<lower=1> D; // max duration 
  int<lower=1> num_coeff; // define because it is used throughout the program
  int<lower=1> lag; // fixed number of lag days (arbitrarily chosen) -- coefficients are a function of duration
  int<lower=1> num_lag_coeff; // define because it is used throughout the program 
  int<lower=1> N; // total number of rows in the dataset
  
  int<lower=1> rows_per_strata[floodzip_id]; // case_control_set*durations[floodzip_id]
  int<lower=0> Y[N]; // Y (cases) is a vector of length N 
  matrix<lower=0, upper = 1>[N,num_coeff] X; // X is a matrix with dimensions N x D*(D+1)/2, takes on values 0 or 1
  matrix<lower=0, upper = 1>[N,num_lag_coeff] A; // A is a matrix with dimensions N x lag, takes on values 0 or 1
  vector<lower=0>[N] offset;
  
  matrix[num_coeff,num_coeff] Sigma_d;
  matrix[num_coeff,num_coeff] Sigma_t;
  matrix[num_lag_coeff,num_lag_coeff] Sigma_d2; 
  matrix[num_lag_coeff,num_lag_coeff] Sigma_l;
  
  // N should be exactly divisible by case_control_set to produce an integer
  int<lower=1> sequence[floodzip_id]; // generate a sequence of numbers to indicate the start of a new strata
  
  vector[num_coeff] mu_beta; // mu_beta is a vector of length D*(D+1)/2
  vector[num_lag_coeff] mu_theta; // mu_theta is a vector of length lag 
}

transformed data {
  vector[N] log_offset = log(offset);
}

// The parameters accepted by the model.
parameters {
  real<lower = 0>little_sigma_beta2; // non-zero 
  real<lower = 0>phi; // non-zero, lengthscale for d -- use this in flooded days 
  real<lower = 0>tau; // non-zero, lengthscale for t
  
  real<lower = 0>little_sigma_theta2; // non-zero
  real<lower = 0>gamma; //non-zero, lengthscale for d -- use this in lags 
  real<lower = 0>eta; // non-zero, lengthscale for l 


  vector[num_coeff] beta; // beta is a vector of length D*(D+1)/2
  vector[num_lag_coeff] theta; // theta is a vector of length lag*D (arbitrary number of days added)
}
  
model { 
  matrix[num_coeff,num_coeff] Kd = exp(-Sigma_d/phi);
  matrix[num_coeff,num_coeff] Kt = exp(-Sigma_t/tau);
  matrix[num_coeff, num_coeff] Sigma_beta = little_sigma_beta2*elt_multiply(Kd,Kt);
  
  matrix[num_coeff, num_coeff] Sigma_beta_reg;
  Sigma_beta_reg = Sigma_beta + diag_matrix(rep_vector(1e-8, num_coeff));
  
  matrix[num_coeff, num_coeff] L_Sigma_beta;
  L_Sigma_beta = cholesky_decompose(Sigma_beta_reg);
  
  matrix[num_lag_coeff,num_lag_coeff] Kd2 = exp(-Sigma_d2/gamma);
  matrix[num_lag_coeff,num_lag_coeff] Kl = exp(-Sigma_l/eta);
  matrix[num_lag_coeff, num_lag_coeff] Sigma_theta = little_sigma_theta2*elt_multiply(Kd2,Kl);
  
  matrix[num_lag_coeff, num_lag_coeff] Sigma_theta_reg;
  Sigma_theta_reg = Sigma_theta + diag_matrix(rep_vector(1e-8, num_lag_coeff));
  
  matrix[num_lag_coeff, num_lag_coeff] L_Sigma_theta;
  L_Sigma_theta = cholesky_decompose(Sigma_theta_reg);
  
  
  { vector[N] log_numer;
  for (obs in 1:N){
    log_numer[obs] = X[obs] * beta + A[obs] * theta + log(offset[obs]); 
  }

  
  for (s in 1:floodzip_id){
    int start = sequence[s];
    int stop = start + rows_per_strata[s]-1;
    
    vector[rows_per_strata[s]] eta_s
      = log_numer[start:stop];
      
    vector[rows_per_strata[s]] pi_s = softmax(eta_s); // sums to 1 by construction 
    
    target += multinomial_lupmf(Y[start:stop]| pi_s);
  }
  }

  // priors
  // variance (little_sigma) should use inv_gamma and precision should use gamma (phi, tau, gamma, eta) - 7/11/24
  little_sigma_beta2 ~ lognormal(log(0.3),0.2); // uniform(0,1); // lognormal(log(0.3), 0.2);
  phi ~ lognormal(log(0.3),0.2); // uniform(0,1); // lognormal(log(0.3), 0.2);
  tau ~ lognormal(log(0.3),0.2); // uniform(0,1) // lognormal(log(0.3), 0.2);
  little_sigma_theta2 ~ lognormal(log(0.3),0.2); // uniform(0,1); // lognormal(log(0.3), 0.2);
  gamma ~ lognormal(log(0.3),0.2); // uniform(0,1); // lognormal(log(0.3), 0.2);
  eta ~ lognormal(log(0.3),0.2); // uniform(0,1) // lognormal(log(0.3), 0.2);
  beta ~ multi_normal_cholesky(mu_beta, L_Sigma_beta); // beta is a vector of length D*(D+1)/2
  theta ~ multi_normal_cholesky(mu_theta, L_Sigma_theta); // theta is a vector of length lag*D 
}
