require(here)
require(purrr)
require(magrittr)
require(dplyr)
require(KFAS)
# require(doParallel)

rm(list=ls())

set.seed(4238)

source(here('source','BSP.R'))
source(here('source','setup.R'))
load(here('output','mortality.Rdata'))

# Data preparation
data_list <- list(ita_man = Y_ita_man/N_ita_man,
                  ita_woman = Y_ita_woman/N_ita_woman)

delta <- 1

collector <- list()
collector <- append(collector, list(kernel = matern_kernel,
                                    delta = delta,
                                    age_knots = age_knots,
                                    rep = rep,
                                    country = names(data_list)))

# Model definition
model_list <- lapply(data_list, 
                     . %>% ngp.model(rates = .,
                                     delta = delta))

# Model fit
fit_list <- lapply(model_list,
                   . %>%
                     ngp.fit(., 
                             rep = 1, 
                             method = 'Nelder-Mead', 
                             parallel = FALSE, 
                             maxcl = 1))

smooth_list <- lapply(fit_list, . %>% 
                        pluck('fit') %>%
                        KFS(., smoothing = c('mean', 'state', 'signal'),
                            simplify = FALSE, 
                            maxiter = 200))

collector <- append(collector, list(warnings = warnings()))

save(list = c('data_list',
              'fit_list',
              'smooth_list',
              'collector'),
     file = here('output', paste('ITA_NGP_delta', delta, '.Rdata', sep = '')))
