######################################
# Helper functions for making plots 


help_weight <- function(x){
  if(sub("^([[:alpha:]]*).*", "\\1", x)=='RE') return(x)
  else return(paste('U',sub(".*?([0-9]+).*", "\\1", x),sep=''))
}


######################################
# KFAS smoothing -> tibble ###########
######################################
# 
# INPUT:
# - out_kfas: output of KFAS::KFS
# - country: country being processed
#
# OUTPUT:
# - tibble of the smoothing distributino
#
smoothing2tibble <- function(out_kfas, country){
  if(country == 'ITA') years <- years[1:(length(years)-1)]
  states_s <- out_kfas$alphahat
  var_states_s <- array_tree(out_kfas$V, margin=3)
  sigma_states_s <- lapply(var_states_s, function(x){sqrt(diag(x))})
  names(sigma_states_s) <- years
  
  sigma_states <- as_tibble(sigma_states_s) %>%
    mutate(state = rep(c('U', 'dU', 'a'), K+1),
           weight = paste('U', rep(0:K, each = 3), sep = '')) %>%
    pivot_longer(-c('state', 'weight'), names_to = 't') %>%
    mutate_at(vars(t), ~ as.numeric(.)) %>%
    mutate(what = 'smoothed')
  
  data_states <- as_tibble(states_s) %>% 
    add_column(what = 'smoothed') %>%
    add_column(t = years) %>%
    pivot_longer(-c(t,what), names_to = 'weight', values_to = 'value') %>%
    mutate(state = sub("^([[:alpha:]]*).*", "\\1", weight)) %>%
    mutate_at(vars(weight), ~ sapply(., help_weight)) %>%
    # Joining the sigma dataset
    inner_join(sigma_states, by = c('state','weight','t','what'), suffix = c("",".sd")) %>%
    mutate_at(vars(weight), ~ factor(., levels= paste('U', rep(0:K), sep = '')))
  
  return(data_states)
}



######################################
# KFAS signal smoothing -> tibble ####
######################################
# 
# INPUT:
# - out_kfas: output of KFAS::KFS
# - country: country being processed
#
# OUTPUT:
# - tibble of the smoothing distributino
#
signal_smoothing2tibble <- function(out_kfas, country){
  if(country == 'ITA') years <- years[1:(length(years)-1)]
  states_s <- out_kfas$alphahat
  d_states_s <- states_s[, paste('dU',0:K, sep = '')]
  d_signal_s <- cbind(d_states_s[,'dU0'], d_states_s[,-1] %*% t(S))
  colnames(d_signal_s) <- ages
  
  
  var_states_s <- array_tree(out_kfas$V, margin=3)
  var_du_splines_s <- lapply(var_states_s, 
                             function(x){
                               x[seq(2+3, nrow(x), by = 3),seq(2+3, ncol(x), by = 3)]})
  var_dsignal_splines_s <- lapply(var_du_splines_s, 
                                  function(x){S %*% x %*% t(S)})
  var_dsignal_0_s <- lapply(var_states_s, 
                            function(x){x[2,2]})
  sigma_dsignal_splines_s <- lapply(var_dsignal_splines_s, function(x){sqrt(diag(x))})
  sigma_dm_s <- map2(var_dsignal_0_s, sigma_dsignal_splines_s, ~ c(sqrt(.x), .y))
  names(sigma_dm_s) <- years
  
  signals_s <- out_kfas$muhat
  var_signals_s <- array_tree(out_kfas$V_mu, margin=3)
  sigma_signals_s <- lapply(var_signals_s, function(x){sqrt(diag(x))})
  names(sigma_signals_s) <- years
  
  sigma_signals <- as_tibble(sigma_signals_s) %>%
    mutate(age = ages) %>%
    pivot_longer(-(age), names_to = 't') %>%
    mutate_at(vars(t), ~ as.numeric(.)) %>%
    mutate(what = 'smoothed',
           state = 'm') %>%
    bind_rows(as_tibble(sigma_dm_s) %>%
                mutate(age = ages) %>%
                pivot_longer(-(age), names_to = 't') %>%
                mutate_at(vars(t), ~ as.numeric(.)) %>%
                mutate(what = 'smoothed',
                       state = 'dm'))
  
  data_signals <- as_tibble(signals_s) %>%
    mutate(what = 'smoothed',
           state = 'm') %>%
    mutate(t = rep(years, 1)) %>%
    pivot_longer(-c(t,what,state), names_to = 'age', values_to = 'value') %>%
    mutate_at(vars(age), as.numeric) %>%
    bind_rows(as_tibble(d_signal_s) %>%
                mutate(what = 'smoothed',
                       state = 'dm') %>%
                mutate(t = rep(years, 1)) %>%
                pivot_longer(-c(t,what,state), names_to = 'age', values_to = 'value') %>%
                mutate_at(vars(age), as.numeric)) %>%
    # Joining the sigma dataset
    inner_join(sigma_signals, by = c('t','what','age','state'), suffix = c("",".sd"))
  return(data_signals)
}

state_comparison_smoothing <- function(fit_lasp, 
                                       year2 = 2020, 
                                       year1 = 2018, 
                                       nsim = 500){
  fit_kfas <- fit_lasp$fit
  sim_states <- simulateSSM(fit_kfas, 
                            type = 'states',
                            filtered = FALSE,
                            conditional = TRUE,
                            nsim = nsim)
  
  index2 <- which(years == year2)
  index1 <- which(years == year1)
  diff <- sim_states[index2,,] - sim_states[index1,,]
  diffU <- diff[paste('U', 0:K, sep = ''),]
  diffdU <- diff[paste('dU', 0:K, sep = ''),]
  
  diff_states <- as_tibble(diffU) %>%
    mutate(spline = paste('S', 0:K, sep = '')) %>%
    pivot_longer(-spline, names_to = "sample") %>%
    mutate(state = 'U') %>%
    bind_rows(as_tibble(diffdU) %>%
                mutate(spline = paste('S', 0:K, sep = '')) %>%
                pivot_longer(-spline, names_to = "sample") %>%
                mutate(state = 'dU'))
  return(diff_states)
}


######################################
# DIFFERENCE LATENT STATES BETWEEN   #
# YEARS2 and                         #
# avg(YEARS1, YEARS1-1, ..., YEARS1-4)
######################################
# 
# INPUT:
# - fit_lasp: model fit
# -  year2
# - year1
#
# OUTPUT:
# - tibble with the difference of the 
#   latent states between year2 and 
#   last 5 years (from year1) average.
#
state_comparison_avg5_smoothing <- function(fit_lasp, 
                                            year2 = 2020, 
                                            year1 = 2019, 
                                            nsim = 500){
  fit_kfas <- fit_lasp$fit
  sim_states <- simulateSSM(fit_kfas, 
                            type = 'states',
                            filtered = FALSE,
                            conditional = TRUE,
                            nsim = nsim)
  
  index2 <- which(years == year2)
  index1 <- which(years %in% c((year1-4):year1))
  avg5 <- apply(sim_states[index1,,], MARGIN = c(2,3), mean)
  diff <- sim_states[index2,,] - avg5
  diffU <- diff[paste('U', 0:K, sep = ''),]
  diffdU <- diff[paste('dU', 0:K, sep = ''),]
  
  diff_states <- as_tibble(diffU) %>%
    mutate(spline = paste('S', 0:K, sep = '')) %>%
    pivot_longer(-spline, names_to = "sample") %>%
    mutate(state = 'U') %>%
    bind_rows(as_tibble(diffdU) %>%
                mutate(spline = paste('S', 0:K, sep = '')) %>%
                pivot_longer(-spline, names_to = "sample") %>%
                mutate(state = 'dU'))
  return(diff_states)
}


######################################
# StMoMo models postprocessing    ####
######################################
# 
# INPUT:
# - item: model forecast 
# - model: label of the model
# - rates: obs data
# - country: which country
#
# OUTPUT:
# - tibble with forecasts 
#
processStMoMo <- function(item, 
                          model, 
                          rates, 
                          country){
  if(country == 'ITA') years <- years[1:(length(years)-1)]
  
  data_obs <- as_tibble(rates) %>%
    mutate(time_pred = years) %>%
    pivot_longer(-time_pred, names_to = 'age', values_to = 'rate.obs') %>%
    mutate_at(vars(age), as.numeric)
  
  item_point_age <- item %>%
    modify(. %>% apply(., MARGIN = c(1,2), mean) %>%
             as_tibble() %>%
             mutate(age=ages) %>%
             pivot_longer(-age, names_to='time_pred', values_to='fit') %>%
             mutate_at(vars(time_pred), as.numeric)) %>%
    map_dfr(.f = ~., .id='time_fit') %>%
    mutate_at(vars(time_fit), . %>% as.numeric(.) %>% `+`(., years[1]+train-2)) %>%   
    mutate(h_ahead = time_pred - time_fit) %>%
    inner_join(data_obs %>% filter(time_pred > years[1]+train-1), by = c('age', 'time_pred')) %>%
    # Building up utilities
    mutate(rmse = (fit-rate.obs)^2,
           mad = abs(fit-rate.obs),
           rmse_log = (log(fit)-log(rate.obs))^2,
           mad_log = abs(log(fit) - log(rate.obs))) %>%
    mutate(model = model)
  
  return(item_point_age)
}

######################################
# HU model postprocessing    #########
######################################
# 
# INPUT:
# - item: model forecast 
# - rates: obs data
# - country: which country
#
# OUTPUT:
# - tibble with forecasts 
#
processHU <- function(item, 
                      rates, 
                      country){
  if(country == 'ITA') years <- years[1:(length(years)-1)]
  
  data_obs <- as_tibble(rates) %>%
    mutate(time_pred = years) %>%
    pivot_longer(-time_pred, names_to = 'age', values_to = 'rate.obs') %>%
    mutate_at(vars(age), as.numeric)
  
  HU_point_age <- item %>%
    modify(. %>% 
             pluck('predicted.values') %>%
             as_tibble() %>%
             mutate(age = ages) %>%
             pivot_longer(-age, names_to = 'time_pred', values_to = 'fit') %>%
             mutate_at(vars(time_pred), as.numeric)) %>%
    map_dfr(.f = ~., .id='time_fit') %>%
    mutate_at(vars(time_fit), . %>% as.numeric(.) %>% `+`(., years[1]+train-2)) %>% 
    inner_join(data_obs, by = c('age', 'time_pred')) %>%
    mutate(h_ahead = time_pred - time_fit) %>%
    # Building up utilities
    mutate(rmse = (fit-rate.obs)^2,
           mad = abs(fit-rate.obs),
           rmse_log = (log(fit)-log(rate.obs))^2,
           mad_log = abs(log(fit) - log(rate.obs))) %>%
    mutate(model = 'HU')
  
  return(HU_point_age)
}


######################################
# CP model postprocessing    #########
######################################
# 
# INPUT:
# - item: model forecast 
# - rates: obs data
# - country: which country
#
# OUTPUT:
# - tibble with forecasts 
#
processCP <- function(item, 
                      rates, 
                      country){
  if(country == 'ITA') years <- years[1:(length(years)-1)]
  
  data_obs <- as_tibble(rates) %>%
    mutate(time_pred = years) %>%
    pivot_longer(-time_pred, names_to = 'age', values_to = 'rate.obs') %>%
    mutate_at(vars(age), as.numeric)
  
  CP_point_age <- item %>%
    modify(. %>% 
             as_tibble(.name_repair = 'unique') %>%
             mutate(age = ages) %>%
             pivot_longer(-age, names_to = 'time_pred', values_to = 'fit') %>%
             mutate_at(vars(time_pred), as.numeric)) %>%
    map_dfr(.f = ~., .id='time_fit') %>%
    drop_na(time_pred) %>%
    mutate_at(vars(time_fit), . %>% as.numeric(.) %>% `+`(., years[1]+train-2)) %>%
    inner_join(data_obs, by = c('age', 'time_pred')) %>%
    mutate(h_ahead = time_pred - time_fit) %>%
    mutate_at(vars(fit), exp) %>%
    # Building up utilities
    mutate(rmse = (fit-rate.obs)^2,
           mad = abs(fit-rate.obs),
           rmse_log = (log(fit)-log(rate.obs))^2,
           mad_log = abs(log(fit) - log(rate.obs))) %>%
    mutate(model = 'CP')
  
  return(CP_point_age)
}


######################################
# BSP model postprocessing    ########
######################################
# 
# INPUT:
# - item: model forecast 
# - rates: obs data
# - country: which country
#
# OUTPUT:
# - tibble with forecasts 
#
processBSP <- function(item, 
                       rates,
                       country){
  if(country == 'ITA') years <- years[1:(length(years)-1)]
  
  data_obs <- as_tibble(rates) %>%
    mutate(time_pred = years) %>%
    pivot_longer(-time_pred, names_to = 'age', values_to = 'rate.obs') %>%
    mutate_at(vars(age), as.numeric)
  
  bsp_point_age <- item %>%
    mutate_at(vars(t), ~ . + years[1]+train-2) %>%
    rename(time_fit = t) %>%
    mutate(time_pred = time_fit + h_ahead) %>% 
    mutate_at(vars(age), as.numeric) %>%
    mutate_at(vars(fit), exp) %>%
    inner_join(data_obs, by = c('age', 'time_pred')) %>%
    # Building up utilities
    mutate(rmse = (fit-rate.obs)^2,
           mad = abs(fit-rate.obs),
           rmse_log = (log(fit)-log(rate.obs))^2,
           mad_log = abs(log(fit) - log(rate.obs))) %>%
    mutate(model = 'BSP')
  
  return(bsp_point_age)
}
