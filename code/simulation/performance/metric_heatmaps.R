library(tidyverse)

simulation <- 
setwd(paste0('/n/holylabs/nethery_lab/Lab/floods-hospitalizations-glm/multinomial_GP/output/simulations/general_simulation', simulation, '/')) #cluster

D <- 
l <- 

#read in metrics 
pct_bias_low25 <- 
pct_bias_xhigh100 <- 
mse_low25 <- 
mse_xhigh100 <- 
coverage_low25 <- 
coverage_xhigh100 <- 

total_coeff <- length(pct_bias_low25) 
prim_coeff <- 1:((D*(D+1))/2)
lag_coeff <- (((D*(D+1))/2)+1):total_coeff

df <- data.frame(d = rep(1:D, each = l),
                 t = rep(1:l, D),
                 pct_bias_low25 = pct_bias_low25[], #fill in [] with prim_coeff or lag_coeff 
                 pct_bias_xhigh100 = pct_bias_xhigh100[],
                 mse_low25 = mse_low25[],
                 mse_xhigh100 = mse_xhigh100[], 
                 coverage_low25 = coverage_low25[],
                 coverage_xhigh100 = coverage_xhigh100[])

df_plot <- data.frame(d = rep(df$d,2), #repeat for the number of settings 
                      t = rep(df$t,2), 
                      pct_bias = c(df$pct_bias_low25, df$pct_bias_xhigh100),
                      mse = c(df$mse_low25, df$mse_xhigh100),
                      coverage = c(df$coverage_low25, df$coverage_xhigh100),
                      type = rep(c("Smooth", "Noisy"), 
                                 each = D*l))

df_plot$type <- factor(df_plot$type, levels = unique(df_plot$type))

pbias_plot <- df_plot %>%
  ggplot(aes(x = t, y = reorder(d, -d), fill = pct_bias)) +
  geom_tile() + 
  facet_wrap(~type) +
  scale_fill_gradient2(low="blue", mid = "white", high="red", midpoint = 0) +
  scale_x_continuous(breaks=seq(1, 14, 1)) +
  labs(x = expression(paste("Day following flood event (",italic("l"),")")), y = expression(paste("Duration (",italic("d"),")")), fill = "Percent bias") +
  theme_bw() +
  theme(legend.position = "bottom")

cov_plot <- df_plot %>%
  ggplot(aes(x = t, y = reorder(d, -d), fill = coverage)) +
  geom_tile() +
  facet_wrap(~type) +
  scale_fill_gradientn(colours = c("blue", "white", "red"),
                       values = scales::rescale(c(0, 0.5, 0.95, 0.975, 1)),
                       limits = c(0,1)) +
  scale_x_continuous(breaks=seq(1, 14, 1)) +
  labs(x = expression(paste("Day following flood event (",italic("l"),")")), y = expression(paste("Duration (",italic("d"),")")), fill = "Coverage probability") +
  theme_bw() +
  theme(legend.position = "bottom")


mse_plot <- df_plot %>%
  ggplot(aes(x = t, y = reorder(d, -d), fill = sqrt(mse))) +
  geom_tile() + 
  facet_wrap(~type) +
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  scale_x_continuous(breaks=seq(1, 14, 1)) +
  labs(x = expression(paste("Day following flood event (",italic("l"),")")), y = expression(paste("Duration (",italic("d"),")")), fill = "RMSE") +
  theme_bw() +
  theme(legend.position = "bottom")
mse_plot

library(ggpubr)
ggarrange(pbias_plot, mse_plot, cov_plot, nrow = 3)
