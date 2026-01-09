library(rstan)
library(tidyverse)

fit <- readRDS(fit_cause.rds) #read in model fit (from stan)

summary_fit <- summary(fit, probs = c(0.025, 0.975))$summary
colnames(model_fit) <- c("mean", "se_mean", "sd", "ll", "ul", "n_eff", "rhat")
model_fit <- as.data.frame(model_fit)

##covariates## 
options(scipen = 999)
zeta_mat_round <- round(model_fit[, c(1,4,5)], 2)  ##subset to coefficients corresponding to covariates## 
write.csv(format(zeta_mat_round, nsmall = 2, scientific = FALSE),
          "zeta_formatted.csv", row.names = FALSE)

##significance and directionality##
model_fit <- model_fit %>% mutate(rr = exp(mean),
                                  rr_ll = exp(ll),
                                  rr_ul = exp(ul))
model_fit <- model_fit %>% mutate(sig = ifelse(ll > 0 | ul < 0, 1, 0),
                                  dir = case_when(mean > 0 & sig == 1 ~ 1, 
                                                  mean < 0 & sig == 1 ~ -1, 
                                                  TRUE ~ 0))


D <- 
df <- data.frame(d = rep(1:D, times = c(1:D)),
                   t = c(1, 1:2 , 1:3, 1:4, 1:5, 1:6, 1:7, 1:8, 1:9, 1:D), 
                   dir = model_fit$dir)

df_plot <- data.frame(d = df$d, 
                      t = df$t, 
                      direction = df$dir)

#make heatmap of coefficient direction
direction_heatmap <- df_plot %>%
  ggplot(aes(x = t, y = reorder(d, -d), fill = as.character(direction))) +
  geom_tile() +
  scale_fill_manual(name = "Direction",
                    values = c("-1" = "blue", "0" = "lightgray", "1" = "red"),
                    labels = c("Negative (Protective)", "Null", "Positive (Harmful)")) +
  scale_x_continuous(breaks=seq(1, 10, 1)) +
  labs(x = expression(paste("Day of flood event (",italic("t"),")")), y = expression(paste("Duration (",italic("d"),")"))) +
  facet_wrap(~ type2) +
  theme_bw() +
  theme(legend.position = "bottom")

