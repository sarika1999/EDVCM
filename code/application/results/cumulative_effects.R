library(tidyverse)

posterior_samples <- read.csv('./multinomial_GP/FFS_final_oct2025/results/cause_2000_2016_covariates_dur10.csv')
dat <- readRDS('./multinomial_GP/FFS_final_oct2025/sampler/cause_2000_2016_spline_covariates_no_lag_dur10_mult_gp_real_data.rds')
covariates <- data.frame(dat[["Z"]])
rows_per_strata <- dat[["rows_per_strata"]]

covariates$strata <- rep(c(1:length(rows_per_strata)), times = rows_per_strata)
durations_by_strata <- covariates %>% group_by(strata) %>% summarise(dur = n()/3)

D <- 10
covariates_split <- vector("list",D)
for (d in 1:D){
  strata_ids <- durations_by_strata %>% filter(dur == d) %>% pull(strata) 
  covariates_split[[d]] <- covariates %>% filter(strata %in% c(strata_ids))
}

# extract posterior samples of interest 
beta_samples <- posterior_samples[, 4:58]   # 55 beta samples
zeta_samples  <- posterior_samples[, 59:76]  # 18 zeta (covariate) coefficients
M <- nrow(beta_samples)

# triangular index for beta 
tri_cols <- lapply(1:D, function(d) {
  start <- d * (d - 1) / 2 + 1
  end   <- d * (d + 1) / 2
  start:end
})

# compute exp(zeta' z_st) for each duration and draw 
# store ??_s exp(zeta' z_st) for each duration and draw

sum_exp_zetaZ <- matrix(0, nrow = M, ncol = D)

for (d in 1:D) { # for each duration 
  Z_d <- as.matrix(covariates_split[[d]] %>% select(-c(19)))   # matrix of covariates for all strata corresponding to floods of duration d
  for (m in 1:M) { #for each row 
    zeta_m <- as.numeric(zeta_samples[m, ]) #get row of covariate coefficient sample 
    linpred <- Z_d %*% zeta_m #multiply coefficient by covariate 
    sum_exp_zetaZ[m, d] <- sum(exp(linpred)) #exponentiate linear predictor and sum over samples (i.e. over all strata)
  }
}

# compute cumulative effect for each duration and row in posterior samples  
delta_by_row <- matrix(0, nrow = M, ncol = D)

for (d in 1:D) { #for each duration 
  cols <- tri_cols[[d]]
  s_t <- sum_exp_zetaZ[, 1:d, drop = FALSE]              # subset 
  betas <- exp(beta_samples[, cols, drop = FALSE])       #exponentiate columns corresponding to a particular duration
  
  num <- rowSums(betas * s_t)                           # sum over t = 1,...,d(exp(beta_{dt})*sum_{i=1}^s \times exp(zeta' z_{st})) #numerator
  den <- rowSums(s_t)                                   # sum over t = 1,..., d #denominator
  delta_by_row[, d] <- num / den
}

# --- Summaries across posterior draws ----
delta_mean <- colMeans(delta_by_row)
delta_ci   <- t(apply(delta_by_row, 2, quantile, probs = c(0.025, 0.975)))

# --- Output results ---
data.frame(
  duration = 1:D,
  mean = delta_mean,
  lower = delta_ci[, 1],
  upper = delta_ci[, 2]
)