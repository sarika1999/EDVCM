# create lags_rand_data.rds 

library(tidyverse)

dat <- readRDS('./data/rand_data.rds')
floodzip_d <- dat %>% select(floodzip_id, d) %>% unique()
lag <- 5  #changeable
case_control_set <- length(unique(dat$control_indicator)) %>% as.double()

last_occurrence <- dat %>% group_by(floodzip_id) %>% filter(t == max(t))
lag_dat <- last_occurrence[rep(seq_len(nrow(last_occurrence)), times = lag), ]

lag_dat <- lag_dat %>% arrange(floodzip_id,zipcode,month,day,control_indicator) #order to match paper 
lag_dat$l <- rep(1:lag, times = nrow(floodzip_d)*case_control_set)
lag_dat$t <- NA
lag_dat$event_exposed <- 0
lag_dat$lag1 <- ifelse(lag_dat$l == 1 & lag_dat$control_indicator == 0, 1, 0)
lag_dat$lag2 <- ifelse(lag_dat$l == 2 & lag_dat$control_indicator == 0, 1, 0)
lag_dat$lag3 <- ifelse(lag_dat$l == 3 & lag_dat$control_indicator == 0, 1, 0)
lag_dat$lag4 <- ifelse(lag_dat$l == 4 & lag_dat$control_indicator == 0, 1, 0)
lag_dat$lag5 <- ifelse(lag_dat$l == 5 & lag_dat$control_indicator == 0, 1, 0)

lag_dat$date = ymd(paste(lag_dat$year, lag_dat$month, lag_dat$day))

lag_dat <- lag_dat %>% mutate(date = case_when(l == 1 ~ date + 1,
                                               l == 2 ~ date + 2,
                                               l == 3 ~ date + 3, 
                                               l == 4 ~ date + 4,
                                               l == 5 ~ date + 5))

lag_dat$year <- year(lag_dat$date)
lag_dat$month <- month(lag_dat$date)
lag_dat$day <- day(lag_dat$date)

#will need to check population counts for floods whose lag days extend into the following year

lag_dat$date <- NULL

dat$l <- NA
dat$lag5 <- dat$lag4 <- dat$lag3 <- dat$lag2 <- dat$lag1 <- 0

lag_dat <- lag_dat[names(dat)]

full_dat <- rbind(dat, lag_dat)
full_dat <- full_dat %>% arrange(floodzip_id,zipcode,month,day,control_indicator) #order to match paper 

saveRDS(full_dat, './data/lags_rand_data.rds')
