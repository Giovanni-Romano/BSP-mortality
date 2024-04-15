setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
require(here)
require(tidyverse)
require(KFAS)
require(parallel)

rm(list=ls())

set.seed(4238)
RNGkind("L'Ecuyer-CMRG")

source(here('source','BSP.R'))
source(here('source','setup.R'))
load(here('output','mortality.Rdata'))

# Data preparation
data_list <- list(ita_man = Y_ita_man/N_ita_man,
                  ita_woman = Y_ita_woman/N_ita_woman,
                  swe_man = Y_swe_man/N_swe_man,
                  swe_woman = Y_swe_woman/N_swe_woman,
                  uk_man = Y_uk_man/N_uk_man,
                  uk_woman = Y_uk_woman/N_uk_woman,
                  us_man = Y_us_man/N_us_man,
                  us_woman = Y_us_woman/N_us_woman)

collector <- list()
collector <- append(collector, list(kernel = matern_kernel,
                                    delta = delta,
                                    age_knots = age_knots,
                                    rep = rep,
                                    country = names(data_list)))

# Model definition
method_lambda <- expand.grid(method.ssm = c("exact", "euler", "3ord"),
                             lambda = c(NA, 1))
model_list <- list()
fit_list <- list()

for (d in seq_along(data_list)) {
  cat("Fitting data", d, "out of", length(data_list), "\n")
  model_list[[d]] <- apply(method_lambda,
                           1,
                           function(m_l)
                             bsp.model(method.ssm = m_l['method.ssm'],
                                       lambda = as.numeric(m_l['lambda']),
                                       rates = data_list[[d]],
                                       delta = delta,
                                       age_knots = age_knots,
                                       kernel = matern_kernel))
  
  names(model_list[[d]]) <- apply(method_lambda %>% 
                                    mutate(lambda = ifelse(is.na(lambda),
                                                           "free",
                                                           "fixed1")),
                                  1,
                                  function(x) paste(x, collapse = '_'))
  
  # Model fit
  fit_list[[d]] <- lapply(model_list[[d]],
                          function(model)
                            bsp.fit(model,
                                    rep = 5,
                                    opt.method = 'Nelder-Mead',
                                    parallel = FALSE,
                                    maxcl = 10))
  
  collector <- append(collector, list(warnings = warnings()))
}

names(model_list) <- names(fit_list) <- names(data_list)


save(list = c('data_list',
              'fit_list',
              'collector'),
     file = here('output', '6_models_fit.Rdata'))


smooth_list <- lapply(fit_list, . %>% 
                        lapply(. %>%
                                 pluck('fit') %>%
                                 KFS(., smoothing = c('mean', 'state', 'signal'),
                                     simplify = FALSE, 
                                     maxiter = 200)))

save(list = c('data_list',
              'fit_list',
              'smooth_list',
              'collector'),
     file = here('output', '6_models_smoothing.Rdata'))