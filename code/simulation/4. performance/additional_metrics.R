library(rstan)
library(tidyverse)

fit <- readRDS(fit.rds) #read in model fit (from stan)

##treedepth##
sp <- get_sampler_params(fit, inc_warmup = FALSE)
max_td <- max(vapply(sp, function(x) max(x[,"treedepth__"]), 0))
hits <- sum(vapply(sp, function(x) sum(x[,"treedepth__"] == max_td), 0))
frac <- hits / sum(vapply(sp, nrow, 0))
c(max_treedepth = max_td, hits = hits, frac = round(frac, 4))

##effective sample size, rhat##
summary_fit <- summary(fit, probs = c(0.025, 0.975))$summary
colnames(model_fit) <- c("mean", "se_mean", "sd", "ll", "ul", "n_eff", "rhat")
model_fit <- as.data.frame(model_fit)