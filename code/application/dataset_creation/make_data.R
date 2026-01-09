rm(list=ls())

# arguments from Rscript
args <- commandArgs(trailingOnly=TRUE)
seedVal = as.numeric(args[1])

# years of analysis
years = c(2000:2016)

library(tidyverse)

# load ccs lookup
code.lookup.merged = read.csv('./floods-hospitalizations-glm/medicare_processing/data_update/CCS_DX.csv')
code.lookup.merged = subset(code.lookup.merged, !(code_chapter%in%c("Complications of pregnancy; childbirth; and the puerperium", 
                                                                    "Congenital anomalies",
                                                                    "Certain conditions originating in the perinatal period",
                                                                    "Symptoms; signs; and ill-defined conditions and factors influencing health status",
                                                                    "Residual codes; unclassified; all E codes [259. and 260.]")))

# make list of broad causes of hospitalization (level 1 names)
causes_groups = unique(as.character(code.lookup.merged$code_chapter))

# process for finding broad causes of death and matching sub causes
causes_group = causes_groups[seedVal]

# directory to load data from
dir.input = paste0('./multinomial_GP/data/real_data/county/expanded_grid_hospitalisations_FFS/',years[1],'_',years[length(years)],'/')

# CCS level 1 input file
dat = readRDS(paste0(dir.input,'medicare_',gsub(" ", "_", causes_group),'_rates_',years[1],'_',years[length(years)],'.rds'))
dat <- dat %>% arrange(floodcty_id, fipscounty, year, month, day, control_indicator)

num_rows <- dat %>% group_by(floodcty_id, control_indicator) %>% summarise(n = n())
num_rows$d <- num_rows$n - 28
dat_10 <- num_rows %>% filter(d <= 10)
dat_subset <- dat %>% filter(floodcty_id %in% dat_10$floodcty_id)

test <- dat_subset %>% group_by(floodcty_id, control_indicator) %>% mutate(r = row_number())

#this is not the best way to do this#
temp <- vector("list", length = length(unique(dat_10$floodcty_id)))
sequence <- seq(1,nrow(dat_10),3)
counter <- 1
for (i in sequence){
  if (i %% 100 == 0){
    print(i)
  }
  temp[[counter]] <- test %>% filter(floodcty_id == dat_10$floodcty_id[i], r <= dat_10$d[i])
  counter <- counter + 1
}


final <- plyr::ldply(temp, data.frame)

#reorder, add indicators

final <- final %>% arrange(floodcty_id, fipscounty, r, control_indicator, day, month, year)
final$t <- final$r

case_control_set <- 3
dur <- final %>% group_by(floodcty_id) %>% summarise(duration = n()/case_control_set)
final <- left_join(final, dur, by = "floodcty_id")

#remove rows with NA for population: counties 46102, 2063
final <- final %>% drop_na()

saveRDS(final, './multinomial_GP/FFS_final_oct2025/medicare_cause_200_2016_no_lag_dur10.rds')

##note: will have to use combine_county_data_cov_lags.R if we want to use covariates or lags in the sampler##
##using covariates is likely to result in a slightly smaller dataset due to missing values##

