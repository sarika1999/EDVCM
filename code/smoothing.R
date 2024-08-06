#number of basis functions 
num_bf <- 3

#1-dimensional
#library(splines)
#beta.range <- seq(-1,1,length.out=1000)
##natural spline for values from 0 to 1 with degrees of freedom equal to the number of basis functions 
#X <- ns(beta.range, df=num_bf)
#coef <- rnorm(num_bf)
#smooth_fn <- X %*% coef
##identify a certain number of coefficients 
#points(beta.range[seq(1,1000,by=100)], smooth_fn[seq(1,1000,by=100)],
#       col="red",pch=16)

# D <- 14
# num_coeff <- D*(D + 1)/2
# Sigma_d <- matrix(NA, nrow = num_coeff, ncol = num_coeff)
# Sigma_t <- matrix(NA, nrow = num_coeff, ncol = num_coeff)
# for (i in 1:num_coeff) {
#   for (j in 1:num_coeff) {
#     d_i <- ceiling((sqrt(1 + 8*i) - 1)/2)
#     d_j <- ceiling((sqrt(1 + 8*j) - 1)/2)
#     t_i <- i - d_i*(d_i-1)/2
#     t_j <- j - d_j*(d_j-1)/2
#     Sigma_d[i,j] <- abs(d_i - d_j)
#     Sigma_t[i,j] <- abs(t_i - t_j)
#   }
# }
# 
# Sigma = exp(-Sigma_d - Sigma_t) #multi-dimensional product without hyperparameters
# 
# mu <- rep(0, num_coeff)
# 
# library(mvtnorm)
# set.seed(2024)
# beta <- as.vector(mvtnorm::rmvnorm(n = 1, mean = mu, sigma = Sigma))

ntimes <- 1:14
d <- rep(1:14, ntimes)
t <- c(1, 1:2, 1:3, 1:4, 1:5, 1:6, 1:7, 1:8, 1:9, 1:10, 1:11, 1:12, 1:13, 1:14)

coord <- matrix(c(d, t), byrow = F, ncol = 2)
colnames(coord) <- c("d", "t")
str(coord)

##https://asbates.rbind.io/2019/03/01/thin-plate-splines/
beta <- rep(1,105) #does not matter, just need something for the dummy model 
dat <- data.frame(coord, beta)

library(lattice)
wireframe(beta ~ d*t, dat) #this looks weird because we do not use expand.grid() i.e. we do not have 1-14 days for each of the 14 durations

#isotropic smooth - bs and m are defaults
#need to understand the inputs 
library(mgcv)

fit <- gam(dat$beta ~ s(dat$d, dat$t, bs = "tp", k = num_bf + 1), fit=FALSE) 

X_2d <- fit$X[,-1] #remove intercept 

#generate magnitude for each basis function from std. normal 
coef <- rnorm(num_bf) #number of basis functions should equal the number of columns in the model matrix (+1 with intercept)

#smooth function is the linear combination of the natural spline values and the magnitudes/coefficients
smooth_fn <- X_2d %*% coef

#plot of smooth function 
plot(smooth_fn)


# install.packages("fields")
# library(fields)
# tps_fields <- Tps(coord, beta, m = num_bf)
# 
# pred_fields <- predict(tps_fields)
# 
# library(tidyverse)
# dat <- dat %>% mutate(fitted_beta = as.vector(pred_fields))
# 
# wire_cloud <- function(x,y,z, point_z,...){
#   panel.wireframe(x,y,z,...)
#   panel.cloud(x,y,point_z,...)
# }
# 
# wireframe(beta ~ d*t, data = dat,
#           panel = wire_cloud, point_z = dat$fitted_beta,
#           pch = 16)



