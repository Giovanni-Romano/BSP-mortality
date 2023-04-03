require(here)
require(tidyverse)
require(MortalityForecast)

rm(list=ls())

set.seed(9382)

load(here('output','mortality_extra.Rdata'))
output_collector <- list()
output_collector <- append(output_collector, list(info=Sys.time()))

# Set up minimum training window, how many step-ahead to predict
train <- 50
h_step <- 10
output_collector <- append(output_collector, list(train = train, h_step = h_step))

# Function to fit  and forecast
fitandforecast <- function(t, train, h_step, obs_rate, years){
  print(t)
  M <- model.HyndmanUllah(data  = t(obs_rate[1:t,]), 
                          x = ages,
                          y = years[1:t])
  P <- predict(M, h = h_step, level = 95)
  return(P)
}

# Function to set up the rolling window: minimum window given by train, then 
# increasing the training set by 1 each time
rolling <- function(cg){
  country <- sub("_.*", "", cg)
  gender <- sub(".*_", "", cg)
  print(paste('Doing', country, gender))
  Y <- eval(parse(text = paste('Y', country, gender, sep = '_')))
  N <- eval(parse(text = paste('N', country, gender, sep = '_')))
  rates <- Y/N
  
  years <- years[1:nrow(Y)]
  
  HUforecast <- lapply(train:length(years),
                       . %>%
                         fitandforecast(t = ., 
                                        h = h_step, 
                                        train = train, 
                                        obs_rate = rates,
                                        years = years))
  
  return(HUforecast)
}

# Set up parameter with country-gender combination
cg <- expand.grid(country = c('fra','dnk','cze'),
                  gender = c('man', 'woman')) %>%
  mutate(cg = paste(country, gender, sep = '_')) %>%
  pull(cg)

# Fit and forecast via lapply
res_forward <- lapply(cg, 
                      rolling)
names(res_forward) <- cg

output_collector <- append(output_collector, list(warnings = warnings()))

save(list = c('output_collector',
              'res_forward'),
     file = here('output','HU_for_extra.Rdata'))
 