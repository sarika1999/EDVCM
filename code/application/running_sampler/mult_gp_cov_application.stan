// April,May 2024
// Multinomial (conditional poisson) with GP prior with mean 'mu_beta' and
// standard deviation 'Sigma_beta' (product of exponential kernels in two dimension).
// with post-flood time (lags) that has a GP prior with mean 'mu_theta' and 
// standard deviation 'Sigma_theta' (exponential kernel in one dimension)

data {
  int<lower=1> floodcty_id; // number of unique flood-zip combinations 
  int<lower=1> case_control_set; // number of cases (2) + controls (1) in a matched set
  vector<lower=1>[floodcty_id] durations; // duration for each unique flood-zip combination 
  int<lower=1> D; // max duration 
  int<lower=1> num_coeff; // define because it is used throughout the program
  int<lower=1> num_conf; // number of covariates/confounders 
  int<lower=1> spline_df; // number of degrees of freedom for spline
  int<lower=1> N; // total number of rows in the dataset
  
  int<lower=1> rows_per_strata[floodcty_id]; // case_control_set*durations[floodcty_id]
  int<lower=0> Y[N]; // Y (cases) is a vector of length N 
  matrix<lower=0, upper = 1>[N,num_coeff] X; // X is a matrix with dimensions N x D*(D+1)/2, takes on values 0 or 1
  matrix[N, num_conf*spline_df] Z; // Z is a matrix with dimensions N x (number of df for spline * number of covariates)
  vector<lower=0>[N] offset;
  
  
  matrix[num_coeff,num_coeff] Sigma_d;
  matrix[num_coeff,num_coeff] Sigma_t;
  
  // N should be exactly divisible by case_control_set to produce an integer
  int<lower=1> sequence[floodcty_id]; // generate a sequence of numbers to indicate the start of a new strata
  
  vector[num_coeff] mu_beta; // mu_beta is a vector of length D*(D+1)/2
}

transformed data {
  vector[N] log_offset = log(offset);   
}

// The parameters accepted by the model.
parameters {
  real<lower = 0> sigma; // non-zero 
  real<lower = 0> phi; // non-zero, lengthscale for d
  real<lower = 0> tau; // non-zero, lengthscale for t
  
  vector[num_coeff] beta; // beta is a vector of length D*(D+1)/2
  vector[num_conf*spline_df] zeta; // zeta is a vector of length covariate/confounders 
}

transformed parameters {
  real<lower=0> little_sigma_beta2 = square(sigma);
}
  
model {
  // matrix[num_coeff, num_coeff] Sigma_beta; // cannot define cov_matrix in model block
  // vector[floodcty_id] log_denom;
  // vector[N] prob;
  // Sigma_beta = exp((-1/phi)*Sigma_d + (-1/tau)*Sigma_t + log(little_sigma_beta2)); // multi-dimensional product 
  
  // matrix[num_coeff, num_coeff] Sigma = little_sigma*exp(-1/phi*Sigma_d + (-1/tau*Sigma_t)); where Sigma_d[i,j] = |d_i - d_j|, Sigma_t[i,j] = |t_i - t_j|
  matrix[num_coeff, num_coeff] Kd  = exp(-Sigma_d / phi);
  matrix[num_coeff, num_coeff] Kt  = exp(-Sigma_t / tau);
  matrix[num_coeff, num_coeff] Sigma_beta = little_sigma_beta2 * elt_multiply(Kd, Kt);

  // jitter for numerical stability
  matrix[num_coeff, num_coeff] Sigma_beta_reg;
  Sigma_beta_reg = Sigma_beta + diag_matrix(rep_vector(1e-8, num_coeff));
  
  matrix[num_coeff, num_coeff] L_Sigma_beta;
  L_Sigma_beta = cholesky_decompose(Sigma_beta_reg);
  
  {vector[N] log_numer;
  for (obs in 1:N){
    log_numer[obs] = X[obs] * beta + Z[obs] * zeta + log(offset[obs]); 
  }
  
  for (s in 1:floodcty_id) {
      int start = sequence[s];
      int stop  = start + rows_per_strata[s] - 1;

      vector[rows_per_strata[s]] eta_s
        = log_numer[start:stop];

      vector[rows_per_strata[s]] pi_s = softmax(eta_s);   // sums to 1 by construction

      target += multinomial_lupmf( Y[start:stop] | pi_s );
    }
  }
  

  // priors
  // variance (little_sigma) should use inv_gamma and precision should use gamma (phi, tau, eta) - 7/11/24
  sigma ~ student_t(3,0,1); // uniform(0,1); // lognormal(log(0.3), 0.2);
  phi ~ lognormal(log(0.3), 0.2); // uniform(0,1)
  tau ~ lognormal(0, 0.6); // uniform(0,1) // lognormal(log(0.3), 0.2);
  zeta ~ normal(0, 100); // independent normal priors on each zeta coefficient (confounder)
  beta ~ multi_normal_cholesky(mu_beta, L_Sigma_beta); // beta is a vector of length D*(D+1)/2
}
