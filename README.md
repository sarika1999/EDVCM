## Multinomial GP

This repository currently contains STAN code and corresponding R scripts for a Multinomial model with a Gaussian Process prior. 
Code is based on the latest version used on the FASSE VDI. 

'mult_gp.stan' : dynamic covariance matrix [hyperparameters are sampled from inv-gamma(5,5)]

'mult_gp_general.stan' : original + considers all rows for a flood-zipcode as a single stratum --> independence of stratum holds 
