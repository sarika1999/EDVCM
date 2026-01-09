setwd('./multinomial_GP/')

indices <- data.frame("simulate" = rep(1:5000, each = 14), "d" = rep(1:14, 5000))
rand_data <- readRDS("./multinomial_GP/data/rand_data.rds")
Y <- readRDS('./multinomial_GP/data/simulations/general_simulation#/rand_Y_general_smooth_med_Xpctnoise_simulation#.rds')

for (seedVal in 1:70000){
  if (seedVal %% 100 == 0){
    print(seedVal)
  }
  
  valY <- indices[seedVal,1]
  valDUR <- indices[seedVal,2]
  
  # for (seedVal in 1:5000){
  # if (seedVal %% 100 == 0){
  #   print(seedVal)
  # }
  # 
  # rand_data$cases <- Y[,seedVal]
  # library(tidyverse)
  # split_data <- split(rand_data, rand_data$d)
  # 
  # colnames(split_data[[1]])[6] <- 'event_exposed1'
  # split_data[[2]] <- split_data[[2]] %>% mutate(event_exposed2 = ifelse(t == 2 & control_indicator == 0, 1, 0),
  #                                   event_exposed1 = ifelse(t == 2, 0, event_exposed)) %>% select(-c(event_exposed)) %>% select(1:5, 12:11, 6:10)
  # 
  # split_data[[3]] <- split_data[[3]] %>% mutate(event_exposed3 = ifelse(t == 3 & control_indicator == 0, 1, 0),
  #                                               event_exposed2 = ifelse(t == 2 & control_indicator == 0, 1, 0),
  #                                               event_exposed1 = ifelse(t == 2 | t == 3, 0, event_exposed)) %>% select(-c(event_exposed)) %>% select(1:5, 13:11, 6:10)
  # 
  # split_data[[4]] <- split_data[[4]] %>% mutate(event_exposed4 = ifelse(t == 4 & control_indicator == 0, 1, 0),
  #                                               event_exposed3 = ifelse(t == 3 & control_indicator == 0, 1, 0),
  #                                               event_exposed2 = ifelse(t == 2 & control_indicator == 0, 1, 0),
  #                                               event_exposed1 = ifelse(t == 2 | t == 3 | t == 4, 0, event_exposed)) %>% select(-c(event_exposed)) %>% select(1:5, 14:11, 6:10)
  # 
  # split_data[[5]] <- split_data[[5]] %>% mutate(event_exposed5 = ifelse(t == 5 & control_indicator == 0, 1, 0),
  #                                               event_exposed4 = ifelse(t == 4 & control_indicator == 0, 1, 0),
  #                                               event_exposed3 = ifelse(t == 3 & control_indicator == 0, 1, 0),
  #                                               event_exposed2 = ifelse(t == 2 & control_indicator == 0, 1, 0),
  #                                               event_exposed1 = ifelse(t == 2 | t == 3 | t == 4 | t == 5, 0, event_exposed)) %>% select(-c(event_exposed)) %>% select(1:5, 15:11, 6:10)
  # 
  # split_data[[6]] <- split_data[[6]] %>% mutate(event_exposed6 = ifelse(t == 6 & control_indicator == 0, 1, 0),
  #                                               event_exposed5 = ifelse(t == 5 & control_indicator == 0, 1, 0),
  #                                               event_exposed4 = ifelse(t == 4 & control_indicator == 0, 1, 0),
  #                                               event_exposed3 = ifelse(t == 3 & control_indicator == 0, 1, 0),
  #                                               event_exposed2 = ifelse(t == 2 & control_indicator == 0, 1, 0),
  #                                               event_exposed1 = ifelse(t == 2 | t == 3 | t == 4 | t == 5 | t == 6, 0, event_exposed)) %>% select(-c(event_exposed)) %>% select(1:5, 16:11, 6:10)
  # 
  # split_data[[7]] <- split_data[[7]] %>% mutate(event_exposed7 = ifelse(t == 7 & control_indicator == 0, 1, 0),
  #                                               event_exposed6 = ifelse(t == 6 & control_indicator == 0, 1, 0),
  #                                               event_exposed5 = ifelse(t == 5 & control_indicator == 0, 1, 0),
  #                                               event_exposed4 = ifelse(t == 4 & control_indicator == 0, 1, 0),
  #                                               event_exposed3 = ifelse(t == 3 & control_indicator == 0, 1, 0),
  #                                               event_exposed2 = ifelse(t == 2 & control_indicator == 0, 1, 0),
  #                                               event_exposed1 = ifelse(t == 2 | t == 3 | t == 4 | t == 5 | t == 6 | t == 7, 0, event_exposed)) %>% select(-c(event_exposed)) %>% select(1:5, 17:11, 6:10)
  # 
  # split_data[[8]] <- split_data[[8]] %>% mutate(event_exposed8 = ifelse(t == 8 & control_indicator == 0, 1, 0),
  #                                               event_exposed7 = ifelse(t == 7 & control_indicator == 0, 1, 0),
  #                                               event_exposed6 = ifelse(t == 6 & control_indicator == 0, 1, 0),
  #                                               event_exposed5 = ifelse(t == 5 & control_indicator == 0, 1, 0),
  #                                               event_exposed4 = ifelse(t == 4 & control_indicator == 0, 1, 0),
  #                                               event_exposed3 = ifelse(t == 3 & control_indicator == 0, 1, 0),
  #                                               event_exposed2 = ifelse(t == 2 & control_indicator == 0, 1, 0),
  #                                               event_exposed1 = ifelse(t == 2 | t == 3 | t == 4 | t == 5 | t == 6 | t == 7 | t == 8, 0, event_exposed)) %>% select(-c(event_exposed)) %>% select(1:5, 18:11, 6:10)
  # 
  # split_data[[9]] <- split_data[[9]] %>% mutate(event_exposed9 = ifelse(t == 9 & control_indicator == 0, 1, 0),
  #                                               event_exposed8 = ifelse(t == 8 & control_indicator == 0, 1, 0),
  #                                               event_exposed7 = ifelse(t == 7 & control_indicator == 0, 1, 0),
  #                                               event_exposed6 = ifelse(t == 6 & control_indicator == 0, 1, 0),
  #                                               event_exposed5 = ifelse(t == 5 & control_indicator == 0, 1, 0),
  #                                               event_exposed4 = ifelse(t == 4 & control_indicator == 0, 1, 0),
  #                                               event_exposed3 = ifelse(t == 3 & control_indicator == 0, 1, 0),
  #                                               event_exposed2 = ifelse(t == 2 & control_indicator == 0, 1, 0),
  #                                               event_exposed1 = ifelse(t == 2 | t == 3 | t == 4 | t == 5 | t == 6 | t == 7 | t == 8 | t == 9, 0, event_exposed)) %>% select(-c(event_exposed)) %>% select(1:5, 19:11, 6:10)
  # 
  # split_data[[10]] <- split_data[[10]] %>% mutate(event_exposed10 = ifelse(t == 10 & control_indicator == 0, 1, 0),
  #                                                 event_exposed9 = ifelse(t == 9 & control_indicator == 0, 1, 0),
  #                                               event_exposed8 = ifelse(t == 8 & control_indicator == 0, 1, 0),
  #                                               event_exposed7 = ifelse(t == 7 & control_indicator == 0, 1, 0),
  #                                               event_exposed6 = ifelse(t == 6 & control_indicator == 0, 1, 0),
  #                                               event_exposed5 = ifelse(t == 5 & control_indicator == 0, 1, 0),
  #                                               event_exposed4 = ifelse(t == 4 & control_indicator == 0, 1, 0),
  #                                               event_exposed3 = ifelse(t == 3 & control_indicator == 0, 1, 0),
  #                                               event_exposed2 = ifelse(t == 2 & control_indicator == 0, 1, 0),
  #                                               event_exposed1 = ifelse(t == 2 | t == 3 | t == 4 | t == 5 | t == 6 | t == 7 | t == 8 | t == 9 | t == 10, 0, event_exposed)) %>% select(-c(event_exposed)) %>% select(1:5, 20:11, 6:10)
  # 
  # split_data[[11]] <- split_data[[11]] %>% mutate(event_exposed11 = ifelse(t == 11 & control_indicator == 0, 1, 0),
  #                                                 event_exposed10 = ifelse(t == 10 & control_indicator == 0, 1, 0),
  #                                                 event_exposed9 = ifelse(t == 9 & control_indicator == 0, 1, 0),
  #                                                 event_exposed8 = ifelse(t == 8 & control_indicator == 0, 1, 0),
  #                                                 event_exposed7 = ifelse(t == 7 & control_indicator == 0, 1, 0),
  #                                                 event_exposed6 = ifelse(t == 6 & control_indicator == 0, 1, 0),
  #                                                 event_exposed5 = ifelse(t == 5 & control_indicator == 0, 1, 0),
  #                                                 event_exposed4 = ifelse(t == 4 & control_indicator == 0, 1, 0),
  #                                                 event_exposed3 = ifelse(t == 3 & control_indicator == 0, 1, 0),
  #                                                 event_exposed2 = ifelse(t == 2 & control_indicator == 0, 1, 0),
  #                                                 event_exposed1 = ifelse(t == 2 | t == 3 | t == 4 | t == 5 | t == 6 | t == 7 | t == 8 | t == 9 | t == 10 | t == 11, 0, event_exposed)) %>% select(-c(event_exposed)) %>% select(1:5, 21:11, 6:10)
  # 
  # split_data[[12]] <- split_data[[12]] %>% mutate(event_exposed12 = ifelse(t == 12 & control_indicator == 0, 1, 0),
  #                                                 event_exposed11 = ifelse(t == 11 & control_indicator == 0, 1, 0),
  #                                                 event_exposed10 = ifelse(t == 10 & control_indicator == 0, 1, 0),
  #                                                 event_exposed9 = ifelse(t == 9 & control_indicator == 0, 1, 0),
  #                                                 event_exposed8 = ifelse(t == 8 & control_indicator == 0, 1, 0),
  #                                                 event_exposed7 = ifelse(t == 7 & control_indicator == 0, 1, 0),
  #                                                 event_exposed6 = ifelse(t == 6 & control_indicator == 0, 1, 0),
  #                                                 event_exposed5 = ifelse(t == 5 & control_indicator == 0, 1, 0),
  #                                                 event_exposed4 = ifelse(t == 4 & control_indicator == 0, 1, 0),
  #                                                 event_exposed3 = ifelse(t == 3 & control_indicator == 0, 1, 0),
  #                                                 event_exposed2 = ifelse(t == 2 & control_indicator == 0, 1, 0),
  #                                                 event_exposed1 = ifelse(t == 2 | t == 3 | t == 4 | t == 5 | t == 6 | t == 7 | t == 8 | t == 9 | t == 10 | t == 11 | t == 12, 0, event_exposed)) %>% select(-c(event_exposed)) %>% select(1:5, 22:11, 6:10)
  # 
  # split_data[[13]] <- split_data[[13]] %>% mutate(event_exposed13 = ifelse(t == 13 & control_indicator == 0, 1, 0),
  #                                                 event_exposed12 = ifelse(t == 12 & control_indicator == 0, 1, 0),
  #                                                 event_exposed11 = ifelse(t == 11 & control_indicator == 0, 1, 0),
  #                                                 event_exposed10 = ifelse(t == 10 & control_indicator == 0, 1, 0),
  #                                                 event_exposed9 = ifelse(t == 9 & control_indicator == 0, 1, 0),
  #                                                 event_exposed8 = ifelse(t == 8 & control_indicator == 0, 1, 0),
  #                                                 event_exposed7 = ifelse(t == 7 & control_indicator == 0, 1, 0),
  #                                                 event_exposed6 = ifelse(t == 6 & control_indicator == 0, 1, 0),
  #                                                 event_exposed5 = ifelse(t == 5 & control_indicator == 0, 1, 0),
  #                                                 event_exposed4 = ifelse(t == 4 & control_indicator == 0, 1, 0),
  #                                                 event_exposed3 = ifelse(t == 3 & control_indicator == 0, 1, 0),
  #                                                 event_exposed2 = ifelse(t == 2 & control_indicator == 0, 1, 0),
  #                                                 event_exposed1 = ifelse(t == 2 | t == 3 | t == 4 | t == 5 | t == 6 | t == 7 | t == 8 | t == 9 | t == 10 | t == 11 | t == 12 | t == 13, 0, event_exposed)) %>% select(-c(event_exposed)) %>% select(1:5, 23:11, 6:10)
  # 
  # split_data[[14]] <- split_data[[14]] %>% mutate(event_exposed14 = ifelse(t == 14 & control_indicator == 0, 1, 0),
  #                                                 event_exposed13 = ifelse(t == 13 & control_indicator == 0, 1, 0),
  #                                                 event_exposed12 = ifelse(t == 12 & control_indicator == 0, 1, 0),
  #                                                 event_exposed11 = ifelse(t == 11 & control_indicator == 0, 1, 0),
  #                                                 event_exposed10 = ifelse(t == 10 & control_indicator == 0, 1, 0),
  #                                                 event_exposed9 = ifelse(t == 9 & control_indicator == 0, 1, 0),
  #                                                 event_exposed8 = ifelse(t == 8 & control_indicator == 0, 1, 0),
  #                                                 event_exposed7 = ifelse(t == 7 & control_indicator == 0, 1, 0),
  #                                                 event_exposed6 = ifelse(t == 6 & control_indicator == 0, 1, 0),
  #                                                 event_exposed5 = ifelse(t == 5 & control_indicator == 0, 1, 0),
  #                                                 event_exposed4 = ifelse(t == 4 & control_indicator == 0, 1, 0),
  #                                                 event_exposed3 = ifelse(t == 3 & control_indicator == 0, 1, 0),
  #                                                 event_exposed2 = ifelse(t == 2 & control_indicator == 0, 1, 0),
  #                                                 event_exposed1 = ifelse(t == 2 | t == 3 | t == 4 | t == 5 | t == 6 | t == 7 | t == 8 | t == 9 | t == 10 | t == 11 | t == 12 | t == 13 | t == 14, 0, event_exposed)) %>% select(-c(event_exposed)) %>% select(1:5, 24:11, 6:10)
  # 
  # library(purrr)
  # 
  # split_data <- map(split_data, ~ arrange(.x, floodzip_id, control_indicator))
  # 
  # for (d in names(split_data)){
  #   saveRDS(split_data[[d]], file = paste0("./multinomial_GP/data/simulations/frequent_comparison_Xpctnoise/rand_data",seedVal, "_dur",d,".rds"))
  # }
  # }
  
  # directory to load data from
  dir.input = paste0("./multinomial_GP/data/simulations/frequent_comparison_Xpctnoise/")
  
  # CCS level 1 input file
  input.file = paste0(dir.input,'rand_data', valY, '_dur',valDUR,'.rds')
  
  # check to see if a file exists for the analysis
  if(file.exists(input.file)){
    dat = readRDS(paste0(input.file))
  }
  
  dat$logpt <- log(dat$population)
  
  library(gnm) ; library(splines) ; library(dlnm)
  
  model_eqn_str <- paste("cases ~", paste(paste0('event_exposed', 1:valDUR), collapse = " + "))
  model_eqn <- as.formula(model_eqn_str)
  model_eqn 
  
  system.time
  ({
    
    mod_dlm_unconstrained = gnm(model_eqn, 
                                data=dat, 
                                offset=logpt, 
                                eliminate=factor(floodzip_id), family=quasipoisson)
    
  })
  
  
  # provide summary for unadjusted model
  t_est_mean_unadj = as.data.frame(mod_dlm_unconstrained$coefficients[1:valDUR])
  t_est_uncertainty_unadj = confint.default(mod_dlm_unconstrained)[c(1:valDUR), ,drop = FALSE]
  
  dat.results_unadj = data.frame(t=c(1:valDUR),
                                 rr=t_est_mean_unadj,rr.ll=t_est_uncertainty_unadj[,1],rr.ul=t_est_uncertainty_unadj[,2])
  
  rownames(dat.results_unadj) = seq(nrow(dat.results_unadj))
  names(dat.results_unadj)[2] = 'rr'
  
  # output directory
  dir.output.model.summary = paste0('./multinomial_GP/output/simulations/freq/med_smooth_spline_Xpctnoise/results/unconstrained_dlm/summary/')
  ifelse(!dir.exists(dir.output.model.summary), dir.create(dir.output.model.summary, recursive=TRUE), FALSE)

  # save model summary 
  write.csv(dat.results_unadj, paste0(dir.output.model.summary,'rand_data', valY, '_dur',valDUR,'_model_summary.csv'))
}
