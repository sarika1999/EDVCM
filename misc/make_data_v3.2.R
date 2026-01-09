library(tidyverse)

dat <- readRDS('./multinomial_GP/data/events_with_matched_controls_nolag_df_v3.rds')
length(unique(dat$floodzip_id)) #46680 flood-zip ids #sanity check 

dat <- dat %>% arrange(floodzip_id,zipcode,month,day,control_indicator) #order to match paper 

events <- dat %>% filter(control_indicator == 0) %>% group_by(floodzip_id) %>% summarise(dur = max(days_since_flood))
table(events$dur) #over half of the flood-zip ids are present in durations 1-6 (24007/46680)

D <- 42

duration <- dat %>% group_by(floodzip_id) %>% summarise(max = max(t))

exposure <- matrix(0, nrow = length(unique(dat$floodzip_id))*3, ncol = D + 2)

for (j in 1:42){
  print(j)
  for (i in 1:length(unique(dat$floodzip_id))){
    exposure[i,j] <- ifelse(duration$max[i] >= j, 1, exposure[i,j])
  }
}

exposure[,43] <- unique(dat$floodzip_id)
exposure[,44] <- 0
exposure[46681:nrow(exposure),1:42] <- 0
exposure[46681:nrow(exposure),43] <- rep(unique(dat$floodzip_id), 2)
exposure[46681:nrow(exposure),44] <- rep(c(1,2), each = length(unique(dat$floodzip_id)))

exposure <- as.data.frame(exposure) 
colnames(exposure) <- c(paste0("day", 1:D), "floodzip_id", "control_indicator") 
exposure$control_indicator <- as.numeric(exposure$control_indicator) 

exposure <- exposure %>% arrange(floodzip_id) 
exposure <- exposure[,c(43,44,1:42)]

#maneuver cases and offset 
outcome <- matrix(0, nrow = nrow(duration)*3, ncol = 42)

cases_by_floodzip <- dat %>% group_by(floodzip_id, control_indicator) %>% summarise(case = list(cases))
cases_by_floodzip$case <- lapply(cases_by_floodzip$case, `length<-`, max(lengths(cases_by_floodzip$case)))

outcome <- cases_by_floodzip %>% unnest_wider(case, names_sep="_")
colnames(outcome) <- c("floodzip_id", "control_indicator", paste0("case_day", 1:D))

offset <- matrix(0, nrow = nrow(duration)*3, ncol = 42)

populations_by_floodzip <- dat %>% group_by(floodzip_id, control_indicator) %>% summarise(pop = list(population))
populations_by_floodzip$pop <- lapply(populations_by_floodzip$pop, `length<-`, max(lengths(populations_by_floodzip$pop)))

offset <- populations_by_floodzip %>% unnest_wider(pop, names_sep="_")
colnames(offset) <- c("floodzip_id", "control_indicator", paste0("population_day", 1:D))

#remove year-month-day, days_since_flood and join everything together 
dat2 <- dat %>% distinct(floodzip_id, zipcode, control_indicator) %>% inner_join(.,exposure) %>% inner_join(.,outcome) %>% inner_join(.,offset)

#reformatted 
saveRDS(dat2, "./multinomial_GP/data/events_with_matched_controls_nolag_df_v3.2.rds")
