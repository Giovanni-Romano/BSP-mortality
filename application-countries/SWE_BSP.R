require(here)
require(tidyverse)
require(KFAS)
require(parallel)

rm(list=ls())

set.seed(3412)
RNGkind("L'Ecuyer-CMRG")

options(mc.cores = 3)

source(here('source','BSP.R'))
source(here('source','setup.R'))
load(here('output','mortality.Rdata'))

data_list <- list(swe_man = Y_swe_man/N_swe_man,
                  swe_woman = Y_swe_woman/N_swe_woman)

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
     file = here('output', 'SWE_fit.Rdata'))
