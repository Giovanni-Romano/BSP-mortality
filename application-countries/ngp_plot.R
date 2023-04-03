require(here)
require(tidyverse)
require(ggh4x)
require(ggrepel)
require(patchwork)
require(KFAS)
require(splines2)
require(latex2exp)
require(ggflags)
require(plotly)
require(htmlwidgets)
theme_set(
  theme_light() +
    theme(
      strip.background = element_rect(color = 'gray', fill = 'white'),
      strip.text.x = element_text(color = 'black'),
      strip.text.y = element_text(color = 'black')
    )
)

rm(list = ls())

set.seed(4523)

source(here('source','helper_fun.R'))
source(here('source','setup.R'))

load(here('output','mortality.Rdata'))
data_list <- list(ita_man = list(rates = Y_ita_man/N_ita_man),
                  ita_woman = list(rates = Y_ita_woman/N_ita_woman))

Z <- length(ages)
K <- Z-1

lambda_est <- NULL
par_est <- tibble(par = c('sigma2_u','sigma2_a','sigma2_e'))
sigma2_e_est <- tibble()
smoothing_tibble <- tibble()
signal_smoothing_tibble <- tibble()
diff_state_20202018 <- tibble()
diff_state_2020avg5 <- tibble()
country <- "ITA"
load(here('output','ITA_NGP_delta1.Rdata'))

lambda_est <- fit_list %>%
  modify(. %>%
           pluck('info') %>%
           pluck('optim') %>%
           pluck('par') %>%
           exp() %>%
           `[`(.,1)) %>%
  imap_dfr(.f = ~ .)
par_est <- fit_list %>%
  modify(. %>% 
           pluck('info') %>%
           pluck('optim') %>%
           pluck('par') %>%
           exp()) %>%
  imap_dfr(.f = ~ .)
smoothing_tibble <- lapply(smooth_list, . %>% 
                             smoothing2tibble(., country = country)) %>%
  bind_rows(.id = 'data') %>%
  mutate(gender = sub(".*_", "", data),
         country = sub("_.*", "", data)) %>%
  mutate(state_gender = paste(state, gender, sep = '_'))

as_tibble(log((Y_ita_man)/(N_ita_man))) %>%
  mutate(t = years[1:(length(years)-1)]) %>%
  pivot_longer(-t, names_to = 'age', values_to = 'value') %>%
  mutate_at(vars(age), as.numeric) %>%
  mutate(gender = "man",
         coutry = "ita") %>%
  bind_rows(as_tibble(log((Y_ita_woman)/(N_ita_woman ))) %>%
              mutate(t = years[1:(length(years)-1)]) %>%
              pivot_longer(-t, names_to = 'age', values_to = 'value') %>%
              mutate_at(vars(age), as.numeric) %>%
              mutate(gender = "woman",
                     coutry = "ita")) -> data_ita


# smoothing_tibble %>%
#   filter(t %in% c(1944, 1960, 1980),
#          state == "U") %>%
#   mutate(age = as.numeric(str_extract(weight, "\\d+$"))) %>%
#   ggplot(aes(x = age, y = value)) +
#   geom_line() + 
#   geom_point(data = data_ita %>%
#                filter(t %in% c(1944, 1960, 1980))) +
#   facet_grid(cols = vars(gender), rows = vars(t))

smoothing_tibble %>%
  mutate(age = as.numeric(str_extract(weight, "\\d+$"))) %>%
  filter(age %in% c(20, 80),
         state == "U") %>%
  ggplot(aes(x = t, y = value)) +
  geom_line() + 
  geom_ribbon(aes(ymin = value - 2*(value.sd),
                  ymax = value + 2*(value.sd)),
              alpha = 0.3) +
  geom_point(data = data_ita %>%
               filter(age %in% c(20, 80))) +
  facet_grid(cols = vars(gender), rows = vars(age), scale = "free_y")










