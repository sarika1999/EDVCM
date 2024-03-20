#simulation results

library(rstan)
test <- readRDS('~/Desktop/multinomial_gp/output/simulations/allfloodzips_dur3_simulation1.rds')

print(test, pars = c("beta", "little_sigma2", "phi", "tau", "lp__"))
pairs(test, pars = c("beta", "little_sigma2", "phi", "tau", "lp__"))

library(bayesplot)
library(tidyverse)

#combine over each chain #only keeps non-warmup iterations
beta_posterior <- as.data.frame(test) %>% 
                  select(-c("little_sigma2", "phi", "tau", "lp__"))

color_scheme_set("brightblue")
uncertainty <- mcmc_intervals(beta_posterior) + 
               scale_y_discrete(labels=c("beta[1]" = expression(beta[1][1]), 
                                         "beta[2]" = expression(beta[2][1]),
                                         "beta[3]" = expression(beta[2][2]), 
                                         "beta[4]" = expression(beta[3][1]),
                                         "beta[5]" = expression(beta[3][2]), 
                                         "beta[6]" = expression(beta[3][3]))) +
              scale_x_continuous(breaks=seq(9,12,by=0.5))
uncertainty

beta_posterior_renamed <- beta_posterior %>% 
                          rename("beta[1][1]" = "beta[1]",
                                 "beta[2][1]" = "beta[2]",
                                 "beta[2][2]" = "beta[3]",
                                 "beta[3][1]" = "beta[4]",
                                 "beta[3][2]" = "beta[5]",
                                 "beta[3][3]" = "beta[6]")

color_scheme_set("purple")
histograms <-mcmc_hist(beta_posterior_renamed, facet_args = list(labeller=ggplot2::label_parsed))
histograms

          