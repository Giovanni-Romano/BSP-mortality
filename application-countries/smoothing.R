require(here)
require(tidyverse)
require(KFAS)

rm(list=ls())

set.seed(9833)

collector <- list()

countries <- c('UK','US','ITA','SWE')

for(cc in countries){
  smooth_list <- list()
  load(here('output', 
            paste(cc, '_fit.Rdata', sep = '')))
  smooth_list <- c(smooth_list,
                   lapply(fit_list, . %>% 
                            pluck('fit') %>%
                            KFS(., smoothing = c('mean', 'state', 'signal'),
                                simplify = FALSE, 
                                maxiter = 200)))
  collector <- append(collector, list(warnings = warnings()))
  
  save(list = c('smooth_list',
                'collector'),
       file = here('output',
                   paste(cc, '_smoothing.Rdata', sep = '')))
}
