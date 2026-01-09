library(tidyverse)

df <- data.frame(d = rep(1:14, times = c(1:14)),
                 t = c(1, 1:2 , 1:3, 1:4, 1:5, 1:6, 1:7, 1:8, 1:9, 1:10, 1:11, 1:12, 1:13, 1:14), 
                 low = rand_true_beta_general_smooth_med_25pctnoise,
                 xhigh = rand_true_beta_general_smooth_med_100pctnoise)

df_plot <- data.frame(d = rep(df$d,2), t = rep(df$t,2), beta = c(df$low, df$xhigh), type = rep(c("Smooth", "Noisy"), each = 105))
df_plot$type <- factor(df_plot$type, levels = unique(df_plot$type))

ground_truth<- df_plot %>%
  ggplot(aes(x = t, y = reorder(d, -d), fill = beta)) +
  geom_tile() + 
  facet_wrap(~type) +
  scale_fill_gradient2(low="blue", mid = "white", high="red", midpoint = 0) +
  scale_x_continuous(breaks=seq(1, 14, 1)) +
  labs(x = expression(paste("Day of flood event (",italic("t"),")")), y = expression(paste("Duration (",italic("d"),")")), fill = "True value") +
  theme_bw() +
  theme(legend.position = "bottom")