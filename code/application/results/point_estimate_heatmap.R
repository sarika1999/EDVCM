library(tidyverse)
library(rstan)

cause_covariates_2000_2016 <- read_csv(".multinomial_GP/FFS_final_oct2025/results/cause_2000_2016_covariates_dur10.csv")

df <- data.frame(d = rep(1:10, times = c(1:10)),
                 t = c(1, 1:2 , 1:3, 1:4, 1:5, 1:6, 1:7, 1:8, 1:9, 1:10), 
                 skin = colMeans(skin_2000_2016_covariates_dur10)[4:58],
                 nerv = colMeans(nerv_2000_2016_covariates_dur10)[4:58],
                 sig_skin = model_fit_skin$sig[4:58],
                 dir_skin = model_fit_skin$dir[4:58],
                 sig_nerv = model_fit_nerv$sig[4:58],
                 dir_nerv = model_fit_nerv$dir[4:58])

df_plot <- data.frame(d = rep(df$d, 2), 
                      t = rep(df$t, 2), 
                      beta = c(df$nerv, df$skin),
                      type = rep(c("A", "B"),each=55),
                      type2 = rep(c("C", "D"), each=55),
                      significance = c(df$sig_nerv, df$sig_skin),
                      direction = c(df$dir_nerv, df$dir_skin))

beta_heatmap <- df_plot %>%
  ggplot(aes(x = t, y = reorder(d, -d), fill = exp(beta))) +
  geom_tile() + 
  scale_fill_gradient2(low="blue", mid = "white", high="red", midpoint = 1.0) +
  scale_x_continuous(breaks=seq(1, 10, 1)) +
  labs(x = expression(paste("Day of flood event (",italic("t"),")")), y = expression(paste("Duration (",italic("d"),")")), fill = "Estimated RR") +
  facet_wrap(~ type) +
  theme_bw() +
  theme(legend.position = "bottom")
