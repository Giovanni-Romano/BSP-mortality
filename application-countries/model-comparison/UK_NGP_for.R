require(here)
require(tidyverse)
require(KFAS)
require(parallel)

rm(list=ls())

set.seed(348576)
RNGkind("L'Ecuyer-CMRG")

options(mc.cores = 3)

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

# Model definition
source(here('source','BSP.R'))
source(here('source','setup.R'))

delta <- 1

output_collector <- append(output_collector, list(delta = delta,
                                                  rep = rep,
                                                  lik = 'gaussian',
                                                  n_for = n_for))

# Loading forecast and forward-rolling functions
source(here('source','BSP_forecast.R'))

cg <- c('uk_man',
        'uk_woman')

res_forward <- lapply(cg, rolling_ngp)

names(res_forward) <- cg

output_collector <- append(output_collector, 
                           list(country_gender = cg,
                                warnings = warnings()))

save(list = c('output_collector',
              'res_forward'),
     file = here('output',paste('UK_NGP_delta', delta, '_for.Rdata', sep = '')))

