library(tidyverse)


# add lags #
dat <- readRDS('./multinomial_GP/FFS_final_oct2025/medicare_cause_2000_2016_no_lag_dur10.rds')
dat$d <- dat$duration
dat$duration <- NULL
floodcty_d <- dat %>% select(floodcty_id, d) %>% unique()
lag <- 5
case_control_set <- length(unique(dat$control_indicator)) %>% as.double()

last_occurrence <- dat %>% group_by(floodcty_id) %>% filter(t == max(t))
lag_dat <- last_occurrence[rep(seq_len(nrow(last_occurrence)), times = lag), ]

lag_dat <- lag_dat %>% arrange(floodcty_id,county,month,day,control_indicator) #order to match paper 
lag_dat$l <- rep(1:lag, times = nrow(floodcty_d)*case_control_set)
lag_dat$t <- NA
lag_dat$event_exposed <- 0
lag_dat$lag1 <- ifelse(lag_dat$l == 1 & lag_dat$control_indicator == 0, 1, 0)
lag_dat$lag2 <- ifelse(lag_dat$l == 2 & lag_dat$control_indicator == 0, 1, 0)
lag_dat$lag3 <- ifelse(lag_dat$l == 3 & lag_dat$control_indicator == 0, 1, 0)
lag_dat$lag4 <- ifelse(lag_dat$l == 4 & lag_dat$control_indicator == 0, 1, 0)
lag_dat$lag5 <- ifelse(lag_dat$l == 5 & lag_dat$control_indicator == 0, 1, 0)
lag_dat$lag6 <- ifelse(lag_dat$l == 6 & lag_dat$control_indicator == 0, 1, 0)
lag_dat$lag7 <- ifelse(lag_dat$l == 7 & lag_dat$control_indicator == 0, 1, 0)

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
full_dat <- full_dat %>% arrange(floodcty_id,county,month,day,control_indicator) #order to match paper 

saveRDS(full_dat, './multinomial_GP/FFS_final_oct2025/medicare_cause_2000_2016_lags_dur10.rds')

library(tidyverse)
library(arrow)

#import gridmet 
path_gridmet <- '/lego/environmental/meteorology__gridmet/county_daily/'

# List the first 16 Parquet files (sorted by name by default)
gridmet <- list.files(path = path_gridmet, pattern = "\\.parquet$", full.names = TRUE)
gridmet <- head(gridmet, 17)

# Read and concatenate the files efficiently
gridmet_dat <- lapply(gridmet, read_parquet)
gridmet_combined <- do.call(rbind, gridmet_dat)

gridmet_combined$year <- year(gridmet_combined$date)
gridmet_combined$month <- month(gridmet_combined$date)
gridmet_combined$day <- day(gridmet_combined$date)
gridmet_combined$date <- NULL
gridmet_combined$county <- as.numeric(gridmet_combined$county)

#import air pollution 
path_ap <- '/lego/environmental/air_pollution__schwartz/county_daily/'

# List the first 16 Parquet files (sorted by name by default)
ap <- list.files(path = path_ap, pattern = "\\.parquet$", full.names = TRUE)

# Read and concatenate the files efficiently
ap_dat <- lapply(ap, read_parquet)
ap_combined <- do.call(rbind, ap_dat)
ap_combined$year <- year(ap_combined$date)
ap_combined$month <- month(ap_combined$date)
ap_combined$day <- day(ap_combined$date)
ap_combined$date <- NULL
ap_combined$county <- as.numeric(ap_combined$county)


confounders <- inner_join(gridmet_combined[,c("county", "tmmx", "rmax", "vs", "year", "month", "day")], ap_combined, by = c("county", "year", "month", "day"))

final_df <- left_join(full_dat, confounders, by = c("county", "year", "month", "day"))

# Step 1: Define the columns to check for NA
columns_to_check <- names(final_df)[15:20]

# Step 2: Find IDs where any NA appears in confounder columns
ids_with_na <- final_df %>%
  filter(if_any(all_of(columns_to_check), is.na)) %>%
  pull(floodcty_id) %>% 
  unique()

# Step 3: Remove all rows with those IDs
final_df <- final_df %>%
  filter(!floodcty_id %in% ids_with_na)

final_df <- final_df %>% arrange(floodcty_id,county,month,day,control_indicator) #order to match paper 
length(unique(final_df$floodcty_id)) #28917 #27282
length(unique(final_df$county)) #2883 #2878

saveRDS(final_df, './multinomial_GP/FFS_final_oct2025/medicare_cause_2000_2016_covariates_lags_dur10.rds')

## make covariates data without lags ##

library(tidyverse)

final_dat <- final_df %>% filter(is.na(l))
final_dat$l <- NULL

saveRDS(final_dat, 'medicare_cause_2000_2016_no_lag_covariates_dur10.rds') 

