require(here)
require(tidyverse)
require(KFAS)
require(doParallel)

rm(list=ls())

set.seed(32534)

output_collector <- list()
output_collector <- append(output_collector, list(info=Sys.time()))

## Parameters
# n_for: last *n_for* years to use to build forecasting model
# train: minimum training size
# h_step: how many steps-ahead to forecast
n_for <- 25  
train <- 58  
h_step <- 10 
output_collector <- append(output_collector, list(train = train, 
                                                  h_step = h_step))

# Mortality data
load(here('output','mortality.Rdata'))

years <- years[1:nrow(Y_swe_man)]

# Model definition
source(here('source','BSP.R'))
source(here('source','setup.R'))

output_collector <- append(output_collector, list(delta = delta,
                                                  age_knots = age_knots,
                                                  rep = rep,
                                                  lik = 'gaussian',
                                                  kernel = matern_kernel,
                                                  n_for = n_for))

# Loading forecast and forward-rolling functions
source(here('source','BSP_forecast.R'))

cg <- c('swe_man',
        'swe_woman')

res_forward <- lapply(cg, 
                      . %>% rolling_uq(., 
                                       n_for = n_for))
names(res_forward) <- cg

output_collector <- append(output_collector, 
                           list(country_gender = cg,
                                warnings = warnings()))

save(list = c('output_collector',
              'res_forward'),
     file = here('output','SWE_uq3_for.Rdata'))
