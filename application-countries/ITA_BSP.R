require(here)
require(tidyverse)
require(KFAS)
require(doParallel)

rm(list=ls())

set.seed(4238)

source(here('source','BSP.R'))
source(here('source','setup.R'))
load(here('output','mortality.Rdata'))

# Data preparation
data_list <- list(ita_man = Y_ita_man/N_ita_man,
                  ita_woman = Y_ita_woman/N_ita_woman)

collector <- list()
collector <- append(collector, list(kernel = matern_kernel,
                                    delta = delta,
                                    age_knots = age_knots,
                                    rep = rep,
                                    country = names(data_list)))

# Model definition
model_list <- lapply(data_list, 
                     . %>% bsp.model(rates = .,
                                     delta = delta,
                                     age_knots = age_knots,
                                     kernel = matern_kernel))

# Model fit
fit_list <- lapply(model_list,
                   . %>%
                     bsp.fit(., 
                             rep = rep, 
                             method = 'Nelder-Mead', 
                             parallel = TRUE, 
                             maxcl = 50))

collector <- append(collector, list(warnings = warnings()))

save(list = c('data_list',
              'fit_list',
              'collector'),
     file = here('output', 'ITA_fit.Rdata'))

