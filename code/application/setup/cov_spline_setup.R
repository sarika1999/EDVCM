medicare_cause_2000_2016_covariates_dur10 <- readRDS("./multinomial_GP/data/real_data/county/FFS/medicare_cause_2000_2016_covariates_dur10.rds")

library(Hmisc)

get_spline <- function(x) {
  rcspline.eval(x, nk = 5, inclx = FALSE)
}

tmmx <- get_spline(medicare_cause_2000_2016_covariates_dur10$tmmx)
rmax <- get_spline(medicare_cause_2000_2016_covariates_dur10$rmax)
vs <- get_spline(medicare_cause_2000_2016_covariates_dur10$vs)
pm25 <- get_spline(medicare_cause_2000_2016_covariates_dur10$pm25)
o3 <- get_spline(medicare_cause_2000_2016_covariates_dur10$o3)
no2 <- get_spline(medicare_cause_2000_2016_covariates_dur10$no2)

conf <- cbind(tmmx, rmax, vs, pm25, o3, no2)
saveRDS(conf, './multinomial_GP/data/real_data/county/FFS_final_oct2025/spline_covariates_dur10.rds')
