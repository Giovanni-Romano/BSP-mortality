require(here)
require(tidyverse)
require(KFAS)
require(doParallel)

rm(list=ls())

set.seed(4535)

source(here('source','BSP.R'))
source(here('source','setup.R'))
load(here('output','mortality.Rdata'))

data_list <- list(us_man = Y_us_man/N_us_man,
                  us_woman = Y_us_woman/N_us_woman)

collector <- list()
collector <- append(collector, list(kernel = matern_kernel,
                                    delta = delta,
                                    age_knots = age_knots,
                                    rep = rep,
                                    country = names(data_list)))

# Model definition
model_list <- lapply(data_list, . %>% 
                       bsp.model(rates = .,
                                 delta = delta,
                                 age_knots = age_knots,
                                 kernel = matern_kernel))

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
     file = here('output', 'US_fit.Rdata'))
