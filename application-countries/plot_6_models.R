setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
require(here)
require(tidyverse)
require(KFAS)
require(parallel)

rm(list=ls())

RNGkind("L'Ecuyer-CMRG")

load(here('output','6_models_smoothing.Rdata'))
source(here('source','helper_fun.R'))


S <- fit_list[[1]][[1]]$info$S
K <- fit_list[[1]][[1]]$info$K

lambda_est <- fit_list %>%
  modify(. %>%
           modify(. %>%
                    pluck('info') %>%
                    pluck('optim') %>%
                    pluck('par') %>%
                    exp() %>%
                    pluck('lambda')))
par_est <- fit_list %>%
  modify(. %>% 
           modify(. %>% 
                    pluck('info') %>%
                    pluck('optim') %>%
                    pluck('par') %>%
                    exp()))


smoothing_tibble <- 
  lapply(smooth_list, . %>% 
           lapply(., 
                  . %>% smoothing2tibble) %>% 
           bind_rows(.id = 'method_lambda') %>%
           mutate(lambda = sub(".*_", "", method_lambda),
                  method.ssm = sub("_.*", "", method_lambda))) %>%
  bind_rows(.id = 'data') %>% 
  mutate(gender = sub(".*_", "", data),
         country = sub("_.*", "", data))

if (!dir.exists(here('output', 'comparison_6_models/'))){
  dir.create(here('output', 'comparison_6_models/'))
}

plots_smoothing <- lapply(names(data_list),
                          function(x)
                            smoothing_tibble %>% 
                            filter(state == "U",
                                   data == x) %>% 
                            ggplot(aes(x = t, y = value, 
                                       lty = lambda,
                                       col = method.ssm)) +
                            geom_line(alpha = 0.6,
                                      linewidth = 0.2) +
                            scale_linetype_manual(name = "Lambda", 
                                                  labels = c("Fixed = 1", "Estimated"),
                                                  values = c("solid", "dashed")) +
                            scale_color_manual(name = "Model", 
                                               values = 2:4,
                                               labels = c("3rd ord. d.e.", "Euler nGP", "Exact nGP")) +
                            facet_wrap( ~ weight, scales = "free_y") + 
                            xlab("Year") +
                            ylab("Smoothing")) 

lapply(plots_smoothing,
       function(P)
         ggsave(filename = here('output', 
                                paste0('comparison_6_models/', 
                                       P %>% 
                                         pluck('data') %>%
                                         pluck('data') %>% 
                                         unique() %>% 
                                         as.character(),
                                       '_smoothing.pdf')),
                plot = P,
                width = 16, height = 9))



# save(list = c('data_list',
#               'fit_list',
#               'collector'),
#      file = here('output', 'ITA_fit_gio.Rdata'))