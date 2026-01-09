library(tidyverse)

#create a random sample of floodzip_id with duration <= 14

dat <- readRDS("./data/events_with_matched_controls_nolag_df_v3.rds")
durations <- dat %>% group_by(floodzip_id) %>% select(d) %>% distinct()

#this constitutes 91.1611% of original data with all 42 durations 
dat_sub <- dat %>% filter(d <= 14)
#nrow(dat_sub) #809,187
durations_sub <- dat_sub %>% group_by(floodzip_id) %>% select(d) %>% distinct()
#proportions(table(durations_sub$d))*100 #identify relative proportions 

floodzip_sample <- sample(durations_sub$floodzip_id, 0.10*nrow(durations_sub), replace = FALSE)
dat_sub_final <- dat_sub_final %>% filter(floodzip_id %in% floodzip_sample)
#durations_sub_final <- dat_sub_final %>% group_by(floodzip_id) %>% select(d) %>% distinct()
#proportions(table(durations_sub_final$d))*100 #check sampled proportions are similar 

saveRDS(dat_sub_final, './data/rand_data.rds')
