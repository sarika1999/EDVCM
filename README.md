## Multinomial GP

This repository currently contains STAN code and corresponding R scripts for a Multinomial model with a Gaussian Process prior. 
Code is based on the latest version used on the FASSE VDI. 

'mult_gp.stan' : constant covariance matrix [hyperparameters are set equal to one]

'mult_gp2.stan' : dynamic covariance matrix [hyperparameters are sampled from inv-gamma(5,5)]
