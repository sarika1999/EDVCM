// February 2023
// Multinomial (conditional poisson) with GP kernel
// with mean 'mu' and standard deviation 'Sigma'.

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
  
  matrix[num_coeff,num_coeff] Sigma_d;
  matrix[num_coeff.num_coeff] Sigma_t;
  
  // N should be exactly divisible by case_control_set to produce an integer
  int<lower=0> sum_Y[N / case_control_set]; // obtain sum of cases for every strata (i.e three rows)
  int<lower=1> sequence[N / case_control_set]; // generate a sequence of numbers to indicate the start of a new strata (i.e. every three rows)
  
  vector<lower=0>[num_coeff] mu; // mu is a vector of length D*(D+1)/2
}

// The parameters accepted by the model.
parameters {
  real eta; //non-zero
  real phi; // non-zero, lengthscale for d
  real tau; // non-zero, lengthscale for t
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
  matrix[num_coeff, num_coeff] Sigma = little_sigma*exp(-1/phi*Sigma_d)*exp(-1/tau*Sigma_t);
  for (strata in 1:num_elements(sequence)){ //no length in STAN 
    for (obs in (sequence[strata]): (sequence[strata] + 2)){
      // j = 1,...,N indexes flood-zipcode-day
      // Y ~ Poisson(lambda_j) --> log(lambda_j) = X_j * beta + log(offset_j)
      numer[obs] = exp(X[obs] * beta) * offset[obs]; 
      // log(pi_j) = X_j * beta + log(offset_j) - log(sum_k(exp(X_k * beta)*offset_k))
      denom[obs] = sum(exp(numer[(sequence[strata]): (sequence[strata] + 2)]));
      // pi_j = exp(X_j * beta)*offset_j / sum_k(exp(X_k * beta)*offset_k)
      prob[obs] = numer[obs]/denom[obs];
      log_lik[obs] = (tgamma(sum_Y[strata])^(1 / case_control_set))*((prob[obs])^(Y[obs]))*inv(tgamma(Y[obs]));
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

