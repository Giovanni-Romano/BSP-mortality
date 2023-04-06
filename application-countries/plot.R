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
data_list <- list(uk_man = list(rates = Y_uk_man/N_uk_man),
                  uk_woman = list(rates = Y_uk_woman/N_uk_woman),
                  us_man = list(rates = Y_us_man/N_us_man),
                  us_woman = list(rates = Y_us_woman/N_us_woman),
                  swe_man = list(rates = Y_swe_man/N_swe_man),
                  swe_woman = list(rates = Y_swe_woman/N_swe_woman),
                  ita_man = list(rates = Y_ita_man/N_ita_man),
                  ita_woman = list(rates = Y_ita_woman/N_ita_woman))
K <- length(age_knots) + 3
Z <- length(ages)

# Splines
knots <- age_knots
S <- mSpline(ages[-1], knots = knots, degree = 2, intercept = TRUE)
K <- ncol(S)
for (k in 1:K){
  S[,k] <- S[,k]/max(S[,k])
}
colnames(S) <- paste('U', 1:K, sep = '')

# Creating data for plots
lambda_est <- NULL
par_est <- tibble(par = c('lambda','sigma2_u','sigma2_a','sigma2_e'))
sigma2_e_est <- tibble()
smoothing_tibble <- tibble()
signal_smoothing_tibble <- tibble()
diff_state_20202018 <- tibble()
diff_state_2020avg5 <- tibble()
for(country in c('UK', 'US', 'ITA', 'SWE')){
  load(here('output',
            paste(country,'_fit.Rdata', sep = '')))
  load(here('output',
            paste(country,'_smoothing.Rdata', sep = '')))
  print(country)
  lambda_est <- bind_cols(lambda_est,
                          fit_list %>%
                            modify(. %>%
                                     pluck('info') %>%
                                     pluck('optim') %>%
                                     pluck('par') %>%
                                     exp() %>%
                                     `[`(.,1)) %>%
                            imap_dfr(.f = ~ .))
  par_est <- bind_cols(par_est,
                       fit_list %>%
                         modify(. %>% 
                                  pluck('info') %>%
                                  pluck('optim') %>%
                                  pluck('par') %>%
                                  exp()) %>%
                         imap_dfr(.f = ~ .))
  smoothing_tibble <- bind_rows(smoothing_tibble,
                                lapply(smooth_list, . %>% 
                                         smoothing2tibble(., country = country)) %>%
                                  bind_rows(.id = 'data') %>%
                                  mutate(gender = sub(".*_", "", data),
                                         country = sub("_.*", "", data)) %>%
                                  mutate(state_gender = paste(state, gender, sep = '_')))
  signal_smoothing_tibble <- bind_rows(signal_smoothing_tibble,
                                       lapply(smooth_list, . %>% 
                                                signal_smoothing2tibble(., country = country)) %>%
                                         bind_rows(.id = 'data') %>%
                                         mutate(gender = sub(".*_", "", data),
                                                country = sub("_.*", "", data)))
  if(country != 'ITA'){
    diff_state_20202018 <- bind_rows(diff_state_20202018,
                                     lapply(fit_list, 
                                            . %>% 
                                              state_comparison_smoothing(.,
                                                                         year2 = 2020,
                                                                         year1 = 2018,
                                                                         nsim = 500)) %>%
                                       bind_rows(.id = 'data') %>%
                                       mutate(gender = sub(".*_", "", data),
                                              country = sub("_.*", "", data)))
    diff_state_2020avg5 <-  bind_rows(diff_state_2020avg5,
                                      lapply(fit_list, 
                                             . %>% 
                                               state_comparison_avg5_smoothing(.,
                                                                               year2 = 2020,
                                                                               year1 = 2019,
                                                                               nsim = 500)) %>%
                                        bind_rows(.id = 'data') %>%
                                        mutate(gender = sub(".*_", "", data),
                                               country = sub("_.*", "", data)))
  }
  
}


sigma2_e_est <- par_est %>% 
  filter(par == 'sigma2_e')

smoothing_tibble_adj <- smoothing_tibble %>%
  filter(state %in% c('U','dU')) %>%
  mutate(lambda_delta = case_when(state == 'U' ~ 1,
                                  data == 'uk_man' ~ lambda_est$uk_man*delta,
                                  data == 'uk_woman' ~ lambda_est$uk_woman*delta,
                                  data == 'us_man' ~ lambda_est$us_man*delta,
                                  data == 'us_woman' ~ lambda_est$us_woman*delta,
                                  data == 'swe_man' ~ lambda_est$swe_man*delta,
                                  data == 'swe_woman' ~ lambda_est$swe_woman*delta,
                                  data == 'ita_man' ~ lambda_est$ita_man*delta,
                                  data == 'ita_woman' ~ lambda_est$ita_woman*delta)) %>%
  mutate_at(vars(value, value.sd), ~ .*lambda_delta) %>%
  mutate_at(vars(gender), ~ ifelse(. == 'man', 'Male', 'Female')) %>%
  mutate_at(vars(country), .f = toupper)

signal_smoothing_tibble_adj <- signal_smoothing_tibble %>%
  mutate(lambda_delta = case_when(state == 'm' ~ 1,
                                  data == 'uk_man' ~ lambda_est$uk_man*delta,
                                  data == 'uk_woman' ~ lambda_est$uk_woman*delta,
                                  data == 'us_man' ~ lambda_est$us_man*delta,
                                  data == 'us_woman' ~ lambda_est$us_woman*delta,
                                  data == 'swe_man' ~ lambda_est$swe_man*delta,
                                  data == 'swe_woman' ~ lambda_est$swe_woman*delta,
                                  data == 'ita_man' ~ lambda_est$ita_man*delta,
                                  data == 'ita_woman' ~ lambda_est$ita_woman*delta)) %>%
  mutate(sigma2_e = case_when(data == 'uk_man' ~ sigma2_e_est$uk_man,
                              data == 'uk_woman' ~ sigma2_e_est$uk_woman,
                              data == 'us_man' ~ sigma2_e_est$us_man,
                              data == 'us_woman' ~ sigma2_e_est$us_woman,
                              data == 'swe_man' ~ sigma2_e_est$swe_man,
                              data == 'swe_woman' ~ sigma2_e_est$swe_woman,
                              data == 'ita_man' ~ sigma2_e_est$ita_man,
                              data == 'ita_woman' ~ sigma2_e_est$ita_woman)) %>%
  mutate_at(vars(value, value.sd), ~ .*lambda_delta) %>%
  mutate_at(vars(gender), ~ ifelse(. == 'man', 'Male', 'Female')) %>%
  mutate_at(vars(country), .f = toupper)

diff_state_20202018_adj <- diff_state_20202018 %>%
  mutate(lambda_delta = case_when(state == 'U' ~ 1,
                                  data == 'uk_man' ~ lambda_est$uk_man*delta,
                                  data == 'uk_woman' ~ lambda_est$uk_woman*delta,
                                  data == 'us_man' ~ lambda_est$us_man*delta,
                                  data == 'us_woman' ~ lambda_est$us_woman*delta,
                                  data == 'swe_man' ~ lambda_est$swe_man*delta,
                                  data == 'swe_woman' ~ lambda_est$swe_woman*delta)) %>%
  mutate_at(vars(value), ~ .*lambda_delta) %>%
  mutate_at(vars(gender), ~ ifelse(. == 'man', 'Male', 'Female')) %>%
  mutate_at(vars(country), .f = toupper)

diff_state_2020avg5_adj <- diff_state_2020avg5 %>%
  mutate(lambda_delta = case_when(state == 'U' ~ 1,
                                  data == 'uk_man' ~ lambda_est$uk_man*delta,
                                  data == 'uk_woman' ~ lambda_est$uk_woman*delta,
                                  data == 'us_man' ~ lambda_est$us_man*delta,
                                  data == 'us_woman' ~ lambda_est$us_woman*delta,
                                  data == 'swe_man' ~ lambda_est$swe_man*delta,
                                  data == 'swe_woman' ~ lambda_est$swe_woman*delta)) %>%
  mutate_at(vars(value), ~ .*lambda_delta) %>%
  mutate_at(vars(gender), ~ ifelse(. == 'man', 'Male', 'Female')) %>%
  mutate_at(vars(country), .f = toupper)

# Intermediate save
save.image(here('output','plot.Rdata'))
# load(here('output','plot.Rdata'))



#########################################
## Splines basis plot ###################
#########################################
ages_max <- c('U0' = 0, 
              apply(S, 2, which.max))
age_classes <- tibble(weight = 'U0',
                      range = 'Newborn')
for(i in 1:K){
  range0 <- which(S[,i] > 0.2)
  age_classes <- bind_rows(age_classes,
                           tibble(weight = paste('U', i, sep = ''),
                                  range = paste(min(range0), max(range0), sep = '~')))
}
age_classes <- age_classes %>%
  mutate(max = ages_max) %>%
  mutate(spline = paste('S',0:K,sep=''))

colors <- c("ITA" = "#abdda4",
            "US" = "#CE3131",
            "UK" = "#2b83ba",
            "SWE" = "#fdae61")

colors_spline <- c("black" = 'grey80',
                   "color" = '#2c7fb8')

colors_spline_names <- ifelse(paste('U',0:K,sep='') %in% c('U3','U8','U15'), 
                              '#2c7fb8', 
                              'black')

Sfull <- bSpline(seq(1, 100, by = 0.1), 
                 knots = knots, 
                 degree = 2, 
                 intercept = TRUE)
for (k in 1:K){
  Sfull[,k] <- Sfull[,k]/max(Sfull[,k])
}
colnames(Sfull) <- paste('U', 1:K, sep = '')

ages_max_Sfull <- seq(1,100,by=0.1)[apply(Sfull,2,which.max)]


UtoG_Tex <- paste('$g_{',20:1,'}(x)$',sep='') %>%
  sapply(X = .,
         FUN = . %>% latex2exp::TeX(input = .), USE.NAMES = FALSE)

cbind('U0' = rep(0, nrow(S)), S) %>%
  rbind(c(1, rep(0, ncol(S))),.) %>%
  as_tibble() %>%
  mutate(age = ages) %>%
  pivot_longer(-age, names_to = 'weight', values_to = 'value') %>%
  # mutate_at(vars(weight), . %>% paste('S', ., sep = '')) %>%
  mutate_at(vars(value), as.numeric) %>%
  mutate(label = ifelse(weight %in% c('U3','U8','U15'), "color", 'black')) %>%
  ggplot(aes(x = age,
             y = factor(weight, 
                        levels = paste('U', K:0, sep = '')))) +
  geom_tile(aes(fill = label, alpha = value), colour = "white", size = 0.4) +
  labs(y = NULL, x = NULL) +
  guides(fill = "none", alpha = "none") +
  scale_x_continuous(breaks = seq(0,100,by=5), limits = c(-0.5,100)) + #, sec.axis = dup_axis()) +
  scale_y_discrete(position = 'left', labels = UtoG_Tex) +
  scale_fill_manual(values = colors_spline) +
  # scale_fill_grey() +
  # scale_fill_gradient(low = 'white', high = 'black') +
  scale_alpha_continuous(range = c(0,1)) +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size = 9),
        plot.margin = margin(t = 0),
        axis.text.y = element_text(colour = rev(colors_spline_names)),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank()) -> plot_spline

as_tibble(Sfull)  %>%
  mutate(age = seq(1, 100, by = 0.1)) %>%
  pivot_longer(-age, names_to = 'weight', values_to = 'value') %>%
  mutate_at(vars(value), as.numeric) %>%
  mutate(label = ifelse(weight %in% c('U3','U8','U15'), "color", 'black')) %>%
  ggplot(aes(x = age, y = value)) +
  geom_line(aes(group = weight, color = label, size = label), alpha = 0.8) +
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = 1),
               size = 0.3,
               color = colors_spline['black']) +
  geom_point(data = tibble(x = c(0,ages_max_Sfull), y = 1), aes(x=x, y=y),
             shape = 5,
             size = 0.9) +
  geom_label_repel(data = tibble(x = c(0,ages_max_Sfull),
                                 y = 1,
                                 label = ages_max),
                   aes(x = x, y = y, label = label),
                   size = 2,
                   segment.size = 0.2,
                   segment.alpha = 0.5,
                   label.size = 0.1,
                   ylim = c(1.05,NA),
                   min.segment.length = 0.1) +
  scale_size_manual(values = c("black" = 0.4, "color" = 0.8)) +
  scale_color_manual(values = colors_spline) +
  scale_x_continuous(limits = c(-0.5,100)) +
  scale_y_continuous(limits = c(0,1.5)) +
  labs(x = NULL, y = NULL) +
  guides(color = "none", alpha = "none", size = "none") +
  geom_segment(aes(x = 0, y = 0, xend = 100, yend = 0), 
               color = "grey90", size = 1) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(b = 0),
        axis.ticks.y = element_blank()) -> plot_spline2

ggsave(plot_spline2 / plot_spline + 
         plot_layout(height = c(1,2)),
       filename = here('output','spline.pdf'),
       width = 8, height = 4)



#########################################
## Smoothing distribution plot ##########
#########################################

smoothing_tibble_adj %>%
  left_join(age_classes, by = "weight") %>%
  filter(weight %in% c('U0','U7','U11'),
         gender == 'Female') %>%
  mutate(label_dem = case_when(range == 'Newborn' ~ 'Infant (age 0)',
                               range == '23~32' ~ 'Young (age 23 ~ 32)',
                               range == '46~64' ~ 'Adult (age 46 ~ 64)')) %>%
  mutate(alpha = ifelse(country == 'ITA', 0.6, 0.3),
         alpha_rib = ifelse(country == 'ITA', 0.4, 0.2)) %>%
  mutate_at(vars(label_dem), . %>% factor(.,
                                          levels = c('Infant (age 0)',
                                                     'Young (age 23 ~ 32)',
                                                     'Adult (age 46 ~ 64)'))) %>%
  ggplot(aes(x = t, y = value)) + 
  geom_segment(data = tibble(x = 1933,
                             y = 0,
                             xend = 2020,
                             yend = 0,
                             state = 'dU'),
               aes(x = x, y = y, xend = xend, yend = yend),
               alpha = 0.5, linetype = 'dashed') +
  geom_line(aes(color = country),
            alpha = 0.6,
            size = 0.6) +
  geom_ribbon(aes(ymin = value-2*value.sd, ymax = value+2*value.sd,
                  fill = country),  alpha = 0.4) +
  geom_flag(data = tibble(x = c(rep(1943,4),
                                rep(2006,4),
                                rep(1943,4)),
                          gender = 'Female',
                          label_dem = factor(c(rep('Infant (age 0)',4),
                                               rep('Young (age 23 ~ 32)',4),
                                               rep('Adult (age 46 ~ 64)',4)),
                                             levels = c('Infant (age 0)',
                                                        'Young (age 23 ~ 32)',
                                                        'Adult (age 46 ~ 64)')),
                          state = 'U',
                          y = c(-5.2,-5.6,-6,-6.4,
                                -2.3,-2.7,-3.1,-3.5,
                                -5.2,-5.6,-6,-6.4),
                          country = rep(c('it','us','gb','se'),3)),
                     aes(x = x, y = y, country = country), size = 2.5) +
  geom_segment(data = tibble(x = c(rep(1948,4),
                                   rep(2011,4),
                                   rep(1948,4)),
                             y = c(-5.2,-5.6,-6,-6.4,
                                   -2.3,-2.7,-3.1,-3.5,
                                   -5.2,-5.6,-6,-6.4),
                             xend = c(rep(1954,4),
                                      rep(2017,4),
                                      rep(1954,4)),
                             yend = c(-5.2,-5.6,-6,-6.4,
                                      -2.3,-2.7,-3.1,-3.5,
                                      -5.2,-5.6,-6,-6.4),
                             country = rep(c('ITA','US','UK','SWE'),3),
                             gender = 'Female',
                             label_dem = factor(c(rep('Infant (age 0)',4),
                                                  rep('Young (age 23 ~ 32)',4),
                                                  rep('Adult (age 46 ~ 64)',4)),
                                                levels = c('Infant (age 0)',
                                                           'Young (age 23 ~ 32)',
                                                           'Adult (age 46 ~ 64)')),
                             state = 'U'),
               aes(x = x, y = y, xend = xend, yend = yend, color = country)) +
  facet_nested(factor(state, levels = c('U', 'dU'),
                      labels = c(U = latex2exp::TeX(r'($\beta$)'),
                                 dU = latex2exp::TeX(r'($\partial\beta$)'))) ~  gender + label_dem,
               scale = 'free',
               labeller = labeller(.rows = label_parsed),
               strip = strip_nested(
                 text_x = elem_list_text(colour = c("black", "black")),
                 background_x = elem_list_rect(fill = c("grey95", "white")),
                 by_layer_x = TRUE)) +
  ggh4x::force_panelsizes(rows = c(2, 1)) + 
  facetted_pos_scales(y = list(scale_y_continuous(n.breaks = 3), 
                               scale_y_continuous(n.breaks = 3, limits = c(-0.5,0.5)))) +
  scale_x_continuous(breaks = c(1940,1980,2020)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(linetype = NULL,y = NULL,x = NULL) +
  theme(legend.position = "none",
        strip.text.y.right = element_text(angle = 0),
        text = element_text(size = 9),
        plot.margin = margin(b = 0),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(linetype = 'dashed')) -> plot_female

smoothing_tibble_adj %>%
  left_join(age_classes, by = "weight") %>%
  filter(weight %in% c('U0','U7','U11'),
         gender == 'Male') %>%
  mutate(label_dem = case_when(range == 'Newborn' ~ 'Infant (age 0)',
                               range == '23~32' ~ 'Young (age 23 ~ 32)',
                               range == '46~64' ~ 'Adult (age 46 ~ 64)')) %>%
  mutate(alpha = ifelse(country == 'ITA', 0.6, 0.3),
         alpha_rib = ifelse(country == 'ITA', 0.2, 0.2)) %>%
  mutate_at(vars(label_dem), . %>% factor(.,
                                          levels = c('Infant (age 0)',
                                                     'Young (age 23 ~ 32)',
                                                     'Adult (age 46 ~ 64)'))) %>%
  ggplot(aes(x = t, y = value)) + 
  geom_segment(data = tibble(x = 1933,
                             y = 0,
                             xend = 2020,
                             yend = 0,
                             state = 'dU'),
               aes(x = x, y = y, xend = xend, yend = yend),
               # alpha = 0.5,
               linetype = 'dashed') +
  geom_line(aes(color = country),
            alpha = 0.6,
            size = 0.6) +
  geom_ribbon(aes(ymin = value-2*value.sd, ymax = value+2*value.sd,
                  fill = country),  alpha = 0.4) +
  geom_flag(data = tibble(x = c(rep(1943,4),
                                rep(2006,4),
                                rep(1943,4)),
                          gender = 'Male',
                          label_dem = factor(c(rep('Infant (age 0)',4),
                                               rep('Young (age 23 ~ 32)',4),
                                               rep('Adult (age 46 ~ 64)',4)),
                                             levels = c('Infant (age 0)',
                                                        'Young (age 23 ~ 32)',
                                                        'Adult (age 46 ~ 64)')),
                          state = 'U',
                          y = c(-4.9,-5.3,-5.7,-6.1,
                                -2.1,-2.5,-2.9,-3.3,
                                -4.9,-5.3,-5.7,-6.1),
                          country = rep(c('it','us','gb','se'),3)),
            aes(x = x, y = y, country = country), size = 2.5) +
  geom_segment(data = tibble(x = c(rep(1948,4),
                                   rep(2011,4),
                                   rep(1948,4)),
                             y = c(-4.9,-5.3,-5.7,-6.1,
                                   -2.1,-2.5,-2.9,-3.3,
                                   -4.9,-5.3,-5.7,-6.1),
                             xend = c(rep(1954,4),
                                      rep(2017,4),
                                      rep(1954,4)),
                             yend = c(-4.9,-5.3,-5.7,-6.1,
                                      -2.1,-2.5,-2.9,-3.3,
                                      -4.9,-5.3,-5.7,-6.1),
                             country = rep(c('ITA','US','UK','SWE'),3),
                             gender = 'Male',
                             label_dem = factor(c(rep('Infant (age 0)',4),
                                                  rep('Young (age 23 ~ 32)',4),
                                                  rep('Adult (age 46 ~ 64)',4)),
                                                levels = c('Infant (age 0)',
                                                           'Young (age 23 ~ 32)',
                                                           'Adult (age 46 ~ 64)')),
                             state = 'U'), aes(x = x, y = y, xend = xend, yend = yend, color = country)) +
  facet_nested(factor(state, levels = c('U', 'dU'),
                      labels = c(U = latex2exp::TeX(r'($\beta$)'),
                                 dU = latex2exp::TeX(r'($\partial\beta$)'))) ~  gender + label_dem,
               scale = 'free',
               labeller = labeller(.rows = label_parsed),
               strip = strip_nested(
                 text_x = elem_list_text(colour = c("black", "black")),
                 background_x = elem_list_rect(fill = c("grey95", "white")),
                 by_layer_x = TRUE)) +
  force_panelsizes(rows = c(2, 1)) + 
  facetted_pos_scales(y = list(scale_y_continuous(n.breaks = 3), 
                               scale_y_continuous(n.breaks = 3, limits = c(-0.5,0.5)))) +
  scale_x_continuous(breaks = c(1940,1980,2020)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(linetype = NULL,y = NULL,x = NULL) +
  theme(legend.position = "none",
        strip.text.y.right = element_text(angle = 0),
        text = element_text(size = 9),
        plot.margin = margin(t = 0),
        axis.title.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(linetype = 'dashed')) -> plot_male


ggsave(plot_female / plot_male,
       filename = here('output','smoothing.pdf'),
       width = 10,
       height = 5)

## NEW VERTICAL
smoothing_tibble_adj %>%
  left_join(age_classes, by = "weight") %>%
  filter(weight == 'U0') %>%
  mutate(label_dem = 'Infant (age 0)') %>%
  mutate(alpha = ifelse(country == 'ITA', 0.6, 0.3),
         alpha_rib = ifelse(country == 'ITA', 0.4, 0.2)) %>%
  ggplot(aes(x = t, y = value)) + 
  geom_segment(data = tibble(x = 1933,
                             y = 0,
                             xend = 2020,
                             yend = 0,
                             state = 'dU'),
               aes(x = x, y = y, xend = xend, yend = yend),
               alpha = 0.5, linetype = 'dashed') +
  geom_line(aes(color = country),
            alpha = 0.9,
            linewidth = 0.5) +
  geom_ribbon(aes(ymin = value-2*value.sd, ymax = value+2*value.sd,
                  fill = country),  alpha = 0.4) +
  annotate("rect", xmin = 1936, xmax = 1946, ymin = -6.35, ymax = -4.65,
           alpha = 1, color = "gray", fill = "white") +
  geom_flag(data = tibble(x = c(rep(1938,4),
                                rep(1938,4)),
                          gender = c(rep('Female', 4),
                                     rep('Male', 4)),
                          label_dem = factor(c(rep('Infant (age 0)',4),
                                               rep('Infant (age 0)',4))),
                          state = 'U',
                          y = rep(c(-4.9,-5.3,-5.7,-6.1),2),
                          country = rep(c('it','us','gb','se'),2)),
            aes(x = x, y = y, country = country), size = 2.5) +
  geom_segment(data = tibble(x = c(rep(1940,4),
                                   rep(1940,4)),
                             y = rep(c(-4.9,-5.3,-5.7,-6.1),2),
                             xend = c(rep(1944,4),
                                      rep(1944,4)),
                             yend = rep(c(-4.9,-5.3,-5.7,-6.1),2),
                             country = rep(c('ITA','US','UK','SWE'),2),
                             gender = c(rep('Female',4),
                                        rep('Male',4)),
                             label_dem = factor(c(rep('Infant (age 0)',8))),
                             state = 'U'),
               aes(x = x, y = y, xend = xend, yend = yend, color = country)) +
  facet_nested(factor(state, levels = c('U', 'dU'),
                      labels = c(U = latex2exp::TeX(r'($\beta$)'),
                                 dU = latex2exp::TeX(r'($\partial\beta$)'))) ~ label_dem + gender,
               scale = 'free',
               labeller = labeller(.rows = label_parsed),
               strip = strip_nested(
                 text_x = elem_list_text(colour = c("black", "black")),
                 background_x = elem_list_rect(fill = c("grey95", "white")),
                 by_layer_x = TRUE)) +
  ggh4x::force_panelsizes(rows = c(2, 1)) + 
  facetted_pos_scales(y = list(scale_y_continuous(n.breaks = 3), 
                               scale_y_continuous(n.breaks = 3, limits = c(-0.5,0.5)))) +
  scale_x_continuous(breaks = c(1940,1980,2020)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(linetype = NULL,y = NULL,x = NULL) +
  theme(legend.position = "none",
        strip.text.y.right = element_text(angle = 0),
        text = element_text(size = 9),
        plot.margin = margin(b = 0),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(linetype = 'dashed')) -> plot_infant

smoothing_tibble_adj %>%
  left_join(age_classes, by = "weight") %>%
  filter(weight == 'U7') %>%
  mutate(label_dem = 'Young (age 23 ~ 32)') %>%
  mutate(alpha = ifelse(country == 'ITA', 0.6, 0.3),
         alpha_rib = ifelse(country == 'ITA', 0.4, 0.2)) %>%
  ggplot(aes(x = t, y = value)) + 
  geom_segment(data = tibble(x = 1933,
                             y = 0,
                             xend = 2020,
                             yend = 0,
                             state = 'dU'),
               aes(x = x, y = y, xend = xend, yend = yend),
               alpha = 0.5, linetype = 'dashed') +
  geom_line(aes(color = country),
            alpha = 0.9,
            linewidth = 0.5) +
  geom_ribbon(aes(ymin = value-2*value.sd, ymax = value+2*value.sd,
                  fill = country),  alpha = 0.4) +
  annotate("rect", xmin = 1936, xmax = 1946, ymin = -6.4, ymax = -5.1,
           alpha = 1, color = "gray", fill = "white") +
  geom_flag(data = tibble(x = c(rep(1938,4),
                                rep(1938,4)),
                          gender = c(rep('Female', 4),
                                     rep('Male', 4)),
                          label_dem = factor(c(rep('Young (age 23 ~ 32)',4),
                                               rep('Young (age 23 ~ 32)',4))),
                          state = 'U',
                          y = rep(c(-5.3,-5.6,-5.9,-6.2),2),
                          country = rep(c('it','us','gb','se'),2)),
            aes(x = x, y = y, country = country), size = 2.5) +
  geom_segment(data = tibble(x = c(rep(1940,4),
                                   rep(1940,4)),
                             y = rep(c(-5.3,-5.6,-5.9,-6.2),2),
                             xend = c(rep(1944,4),
                                      rep(1944,4)),
                             yend = rep(c(-5.3,-5.6,-5.9,-6.2),2),
                             country = rep(c('ITA','US','UK','SWE'),2),
                             gender = c(rep('Female',4),
                                        rep('Male',4)),
                             label_dem = factor(c(rep('Young (age 23 ~ 32)',8))),
                             state = 'U'),
               aes(x = x, y = y, xend = xend, yend = yend, color = country)) +
  facet_nested(factor(state, levels = c('U', 'dU'),
                      labels = c(U = latex2exp::TeX(r'($\beta$)'),
                                 dU = latex2exp::TeX(r'($\partial\beta$)'))) ~ label_dem + gender,
               scale = 'free',
               labeller = labeller(.rows = label_parsed),
               strip = strip_nested(
                 text_x = elem_list_text(colour = c("black", "black")),
                 background_x = elem_list_rect(fill = c("grey95", "white")),
                 by_layer_x = TRUE)) +
  ggh4x::force_panelsizes(rows = c(2, 1)) + 
  facetted_pos_scales(y = list(scale_y_continuous(n.breaks = 3), 
                               scale_y_continuous(n.breaks = 3, limits = c(-0.5,0.5)))) +
  scale_x_continuous(breaks = c(1940,1980,2020)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(linetype = NULL,y = NULL,x = NULL) +
  theme(legend.position = "none",
        strip.text.y.right = element_text(angle = 0),
        text = element_text(size = 9),
        plot.margin = margin(b = 0),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(linetype = 'dashed')) -> plot_23_32


smoothing_tibble_adj %>%
  left_join(age_classes, by = "weight") %>%
  filter(weight == 'U11') %>%
  mutate(label_dem = 'Adult (age 46 ~ 64)') %>%
  mutate(alpha = ifelse(country == 'ITA', 0.6, 0.3),
         alpha_rib = ifelse(country == 'ITA', 0.4, 0.2)) %>%
  ggplot(aes(x = t, y = value)) + 
  geom_segment(data = tibble(x = 1933,
                             y = 0,
                             xend = 2020,
                             yend = 0,
                             state = 'dU'),
               aes(x = x, y = y, xend = xend, yend = yend),
               alpha = 0.5, linetype = 'dashed') +
  geom_line(aes(color = country),
            alpha = 0.9,
            linewidth = 0.5) +
  geom_ribbon(aes(ymin = value-2*value.sd, ymax = value+2*value.sd,
                  fill = country),  alpha = 0.4) +
  annotate("rect", xmin = 1936, xmax = 1946, ymin = -4.65, ymax = -4,
           alpha = 1, color = "gray", fill = "white") +
  geom_flag(data = tibble(x = c(rep(1938,4),
                                rep(1938,4)),
                          gender = c(rep('Female', 4),
                                     rep('Male', 4)),
                          label_dem = factor(c(rep('Adult (age 46 ~ 64)',4),
                                               rep('Adult (age 46 ~ 64)',4))),
                          state = 'U',
                          y = rep(c(-4.1,-4.25,-4.4,-4.55),2),
                          country = rep(c('it','us','gb','se'),2)),
            aes(x = x, y = y, country = country), size = 2.5) +
  geom_segment(data = tibble(x = c(rep(1940,4),
                                   rep(1940,4)),
                             y = rep(c(-4.1,-4.25,-4.4,-4.55),2),
                             xend = c(rep(1944,4),
                                      rep(1944,4)),
                             yend = rep(c(-4.1,-4.25,-4.4,-4.55),2),
                             country = rep(c('ITA','US','UK','SWE'),2),
                             gender = c(rep('Female',4),
                                        rep('Male',4)),
                             label_dem = factor(c(rep('Adult (age 46 ~ 64)',8))),
                             state = 'U'),
               aes(x = x, y = y, xend = xend, yend = yend, color = country)) +
  facet_nested(factor(state, levels = c('U', 'dU'),
                      labels = c(U = latex2exp::TeX(r'($\beta$)'),
                                 dU = latex2exp::TeX(r'($\partial\beta$)'))) ~ label_dem + gender,
               scale = 'free',
               labeller = labeller(.rows = label_parsed),
               strip = strip_nested(
                 text_x = elem_list_text(colour = c("black", "black")),
                 background_x = elem_list_rect(fill = c("grey95", "white")),
                 by_layer_x = TRUE)) +
  ggh4x::force_panelsizes(rows = c(2, 1)) + 
  facetted_pos_scales(y = list(scale_y_continuous(n.breaks = 3), 
                               scale_y_continuous(n.breaks = 3, limits = c(-0.5,0.5)))) +
  scale_x_continuous(breaks = c(1940,1980,2020)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(linetype = NULL,y = NULL,x = NULL) +
  theme(legend.position = "none",
        strip.text.y.right = element_text(angle = 0),
        text = element_text(size = 9),
        plot.margin = margin(b = 0),
        axis.title.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(linetype = 'dashed')) -> plot_46_64

plot_infant / plot_23_32 / plot_46_64 -> plot_coeff

ggsave(plot_coeff,
       filename = here('output','smoothing.pdf'),
       width = 10,
       height = 8)

#########################################
## Difference 2020 and average last 5   #
## years (Covid plot)                   #
#########################################

diff_state_2020avg5_adj %>%
  group_by(spline, gender, country, state) %>%
  summarise(value.m = mean(value),
            value.sd = sd(value)) %>%
  left_join(age_classes, by = "spline") %>%
  group_by(gender) %>%
  mutate(max = max(value.m + 2*value.sd) + 0.04,
         min = min(value.m - 2*value.sd) - 0.04) %>%
  pivot_longer(cols = c(max,min),
               names_to = "stripe.what", 
               values_to = "stripe.limit") %>%
  mutate(stripe = ifelse(range %in% age_classes$range[seq(1,K+1,by=2)],
                         0.05,0)) %>%
  group_by(gender, country, range) %>%
  filter(gender == 'Male') %>%
  ggplot(aes(x = factor(range, levels=age_classes$range),
             y = value.m,
             group = country)) +
  geom_col(aes(y = stripe.limit, alpha = stripe), 
           position = 'identity', fill = 'grey') +
  geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.5) +
  geom_tile(data = . %>%
              group_by(gender, country, range, state) %>%
              slice(1),
            aes(y = value.m, 
                height = 4*value.sd,
                fill = country),
            width = 0.9,
            position = position_dodge2(1),
            alpha = 0.4) +
  ggflags::geom_flag(data = . %>%
                       filter(state == 'U') %>%
                       group_by(gender, country, range, state) %>%
                       slice(1),
                     aes(country = case_when(country == 'US' ~ 'us',
                                             country == 'UK' ~ 'gb',
                                             country == 'SWE' ~ 'se'),
                         size = size),
                     size = 3.8,
                     position = position_dodge2(0.9)) +
  ggflags::geom_flag(data = . %>%
                       filter(state == 'dU') %>%
                       group_by(gender, country, range, state) %>%
                       slice(1),
                     aes(country = case_when(country == 'US' ~ 'us',
                                             country == 'UK' ~ 'gb',
                                             country == 'SWE' ~ 'se'),
                         size = size),
                     size = 2.5,
                     position = position_dodge2(0.9)) +
  facet_grid(rows = vars(factor(state, 
                                levels = c('U','dU'),
                                labels = c(U = latex2exp::TeX(r'($\beta)'),
                                           dU = latex2exp::TeX(r'($\partial\beta)')))),
             labeller = labeller(.rows = label_parsed),
             cols = vars(gender),
             scale = 'free_y') +
  ggh4x::force_panelsizes(rows = c(2, 1)) + 
  scale_alpha_continuous(range = c(0,0.05)) + 
  scale_y_continuous(expand = c(0,0), n.breaks = 3) +
  scale_fill_manual(values = colors[c('US','SWE','UK')]) +
  scale_color_manual(values = colors[c('US','SWE','UK')]) +
  guides(alpha = "none", color = "none", fill = "none",
         size = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        plot.margin = margin(t = 0),
        strip.text.y.right = element_text(angle = 0),
        axis.title = element_blank()) -> plot_covid_male

diff_state_2020avg5_adj %>%
  group_by(spline, gender, country, state) %>%
  summarise(value.m = mean(value),
            value.sd = sd(value)) %>%
  left_join(age_classes, by = "spline") %>%
  group_by(gender) %>%
  mutate(max = max(value.m + 2*value.sd) + 0.04,
         min = min(value.m - 2*value.sd) - 0.04) %>%
  pivot_longer(cols = c(max,min),
               names_to = "stripe.what", 
               values_to = "stripe.limit") %>%
  mutate(stripe = ifelse(range %in% age_classes$range[seq(1,K+1,by=2)],
                         0.05,0)) %>%
  group_by(gender, country, range) %>%
  filter(gender == 'Female') %>%
  ggplot(aes(x = factor(range, levels=age_classes$range),
             y = value.m,
             group = country)) +
  geom_col(aes(y = stripe.limit, alpha = stripe), 
           position = 'identity', fill = 'grey') +
  geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.5) +
  geom_tile(data = . %>%
              group_by(gender, country, range, state) %>%
              slice(1),
            aes(y = value.m, 
                height = 4*value.sd,
                fill = country),
            width = 0.9,
            position = position_dodge2(1),
            alpha = 0.4) +
  ggflags::geom_flag(data = . %>%
                       filter(state == 'U') %>%
                       group_by(gender, country, range, state) %>%
                       slice(1),
                     aes(country = case_when(country == 'US' ~ 'us',
                                             country == 'UK' ~ 'gb',
                                             country == 'SWE' ~ 'se'),
                         size = size), 
                     size = 3.8,
                     position = position_dodge2(0.9)) +
  ggflags::geom_flag(data = . %>%
                       filter(state == 'dU') %>%
                       group_by(gender, country, range, state) %>%
                       slice(1),
                     aes(country = case_when(country == 'US' ~ 'us',
                                             country == 'UK' ~ 'gb',
                                             country == 'SWE' ~ 'se'),
                         size = size), 
                     size = 2.5,
                     position = position_dodge2(0.9)) +
  facet_grid(rows = vars(factor(state, 
                                levels = c('U','dU'),
                                labels = c(U = latex2exp::TeX(r'($\beta)'),
                                           dU = latex2exp::TeX(r'($\partial\beta)')))),
             labeller = labeller(.rows = label_parsed),
             cols = vars(gender),
             scale = 'free_y') +
  ggh4x::force_panelsizes(rows = c(2, 1)) + 
  scale_alpha_continuous(range = c(0,0.05)) + 
  scale_y_continuous(expand = c(0,0), n.breaks = 3) +
  scale_fill_manual(values = colors[c('US','SWE','UK')]) +
  scale_color_manual(values = colors[c('US','SWE','UK')]) +
  guides(alpha = "none", color = "none", fill = "none") +
  theme(panel.grid = element_blank(),
        plot.margin = margin(b = 0),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.y.right = element_text(angle = 0),
        axis.title = element_blank()) -> plot_covid_female

plot_covid_female / plot_covid_male -> plot_covid
ggsave(plot_covid,
       filename = here('output','2020_avg5.pdf'),
       width = 10,
       height = 7)

## NEW FLAGS VERSION
diff_state_2020avg5_adj %>%
  filter(state == "U") %>%
  group_by(spline, gender, country, state) %>%
  summarise(value.m = mean(value),
            value.sd = sd(value)) %>%
  left_join(age_classes, by = "spline") %>%
  group_by(gender) %>%
  mutate(max = max(value.m + 2*value.sd) + 0.04,
         min = min(value.m - 2*value.sd) - 0.04) %>%
  pivot_longer(cols = c(max,min),
               names_to = "stripe.what", 
               values_to = "stripe.limit") %>%
  mutate(stripe = ifelse(range %in% age_classes$range[seq(1,K+1,by=2)],
                         0.05,0)) %>%
  group_by(gender, country, range) %>%
  filter(gender == 'Male') %>%
  ggplot(aes(x = factor(range, levels=age_classes$range),
             y = value.m,
             group = country)) +
  geom_col(aes(y = stripe.limit, alpha = stripe), 
           position = 'identity', fill = 'grey') +
  geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.5) +
  geom_tile(data = . %>%
              group_by(gender, country, range, state) %>%
              slice(1),
            aes(y = value.m, 
                height = 4*value.sd,
                fill = country),
            width = 0.9,
            position = position_dodge2(1),
            alpha = 0.4) +
  geom_point(data = . %>%
               filter(state == 'U') %>%
               group_by(gender, country, range, state) %>%
               slice(1),
             aes(color = country),
             size = 3.8,
             position = position_dodge2(0.9),
             alpha = 0.8) +
  geom_line(data = . %>%
              filter(state == 'U') %>%
              group_by(gender, country, range, state) %>%
              slice(1),
            aes(color = country),
            position = position_dodge2(0.9),
            alpha = 0.7) +
  annotate("rect", xmin = 15.5, xmax = 17.5, ymin = -0.28, ymax = -0.16,
           alpha = 1, color = "gray", fill = "white") +
  geom_flag(data = tibble(x = rep(16,3),
                          gender = 'Male',
                          y = c(-0.19,-0.22,-0.25),
                          country = c('us','gb','se')),
            aes(x = x, y = y, country = country), 
            size = 4) +
  geom_segment(data = tibble(x = rep(16.4,3),
                             y = c(-0.19,-0.22,-0.25),
                             xend = rep(17,3),
                             yend = c(-0.19,-0.22,-0.25),
                             country = c('US','UK','SWE'),
                             gender = 'Male'),
               aes(x = x, y = y, 
                   xend = xend, yend = yend, 
                   color = country),
               linewidth = 1.8,
               alpha = 0.8) +
  facet_grid(cols = vars(gender),
             scale = 'free_y') +
  ggh4x::force_panelsizes(rows = c(2, 1)) + 
  scale_alpha_continuous(range = c(0,0.05)) + 
  scale_y_continuous(expand = c(0,0), n.breaks = 3) +
  scale_fill_manual(values = colors[c('US','SWE','UK')]) +
  scale_color_manual(values = colors[c('US','SWE','UK')]) +
  guides(alpha = "none", color = "none", fill = "none",
         size = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        plot.margin = margin(t = 0),
        strip.text.y.right = element_text(angle = 0),
        axis.title = element_blank()) -> plot_covid_male


diff_state_2020avg5_adj %>%
  filter(state == "U") %>%
  group_by(spline, gender, country, state) %>%
  summarise(value.m = mean(value),
            value.sd = sd(value)) %>%
  left_join(age_classes, by = "spline") %>%
  group_by(gender) %>%
  mutate(max = max(value.m + 2*value.sd) + 0.04,
         min = min(value.m - 2*value.sd) - 0.04) %>%
  pivot_longer(cols = c(max,min),
               names_to = "stripe.what", 
               values_to = "stripe.limit") %>%
  mutate(stripe = ifelse(range %in% age_classes$range[seq(1,K+1,by=2)],
                         0.05,0)) %>%
  group_by(gender, country, range) %>%
  filter(gender == 'Female') %>%
  ggplot(aes(x = factor(range, levels=age_classes$range),
             y = value.m,
             group = country)) +
  geom_col(aes(y = stripe.limit, alpha = stripe), 
           position = 'identity', fill = 'grey') +
  geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.5) +
  geom_tile(data = . %>%
              group_by(gender, country, range, state) %>%
              slice(1),
            aes(y = value.m, 
                height = 4*value.sd,
                fill = country),
            width = 0.9,
            position = position_dodge2(1),
            alpha = 0.4) +
  geom_point(data = . %>%
               filter(state == 'U') %>%
               group_by(gender, country, range, state) %>%
               slice(1),
             aes(color = country),
             size = 3.8,
             position = position_dodge2(0.9),
             alpha = 0.8) +
  geom_line(data = . %>%
              filter(state == 'U') %>%
              group_by(gender, country, range, state) %>%
              slice(1),
            aes(color = country),
            position = position_dodge2(0.9),
            alpha = 0.7) +
  annotate("rect", xmin = 15.5, xmax = 17.5, ymin = -0.31, ymax = -0.19,
           alpha = 1, color = "gray", fill = "white") +
  geom_flag(data = tibble(x = rep(16,3),
                          gender = 'Female',
                          y = c(-0.22,-0.25,-0.28),
                          country = c('us','gb','se')),
            aes(x = x, y = y, country = country), 
            size = 4) +
  geom_segment(data = tibble(x = rep(16.4,3),
                             y = c(-0.22,-0.25,-0.28),
                             xend = rep(17,3),
                             yend = c(-0.22,-0.25,-0.28),
                             country = c('US','UK','SWE'),
                             gender = 'Female'),
               aes(x = x, y = y, 
                   xend = xend, yend = yend, 
                   color = country),
               linewidth = 1.8,
               alpha = 0.8) +
  facet_grid(cols = vars(gender),
             scale = 'free_y') +
  ggh4x::force_panelsizes(rows = c(2, 1)) + 
  scale_alpha_continuous(range = c(0,0.05)) + 
  scale_y_continuous(expand = c(0,0), n.breaks = 3) +
  scale_fill_manual(values = colors[c('US','SWE','UK')]) +
  scale_color_manual(values = colors[c('US','SWE','UK')]) +
  guides(alpha = "none", color = "none", fill = "none") +
  theme(panel.grid = element_blank(),
        plot.margin = margin(b = 0),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.y.right = element_text(angle = 0),
        axis.title = element_blank()) -> plot_covid_female

plot_covid_female / plot_covid_male -> plot_covid
ggsave(plot_covid,
       filename = here('output','2020_avg5.pdf'),
       width = 10,
       height = 7)

#########################################
## 3-d Data plot      ###################
#########################################

as_tibble(log((Y_ita_man + Y_ita_woman)/(N_ita_man + N_ita_woman))) %>%
  mutate(year = years[1:(length(years)-1)]) %>%
  pivot_longer(-year, names_to = 'age', values_to = 'value') %>%
  mutate_at(vars(age), as.numeric) %>%
  pivot_wider(names_from = age, values_from = value) %>%
  select(-year) %>%
  as.matrix() %>% 
  t() -> data3d_ita
colnames(data3d_ita) <- years[1:(length(years)-1)]

as_tibble(log((Y_us_man + Y_us_woman)/(N_us_man + N_us_woman))) %>%
  mutate(year = years) %>%
  pivot_longer(-year, names_to = 'age', values_to = 'value') %>%
  mutate_at(vars(age), as.numeric) %>%
  pivot_wider(names_from = age, values_from = value) %>%
  select(-year) %>%
  as.matrix() %>% 
  t() -> data3d_us
colnames(data3d_us) <- years


plot_ly(z=~data3d_ita,
        x = years[1:(length(years)-1)]) %>%
  add_surface(alpha = 0.8,
              colors = "Spectral",
              reversescale = TRUE) %>%
  add_trace(x = rep(1970, 101), 
            y = 0:100, 
            z = data3d_ita[,'1970'], 
            mode = 'lines',
            type = 'scatter3d',
            line = list(color = '#037A73',
                        width = 6)) %>%
  add_trace(x = rep(2005, 101), 
            y = 0:100, 
            z = data3d_ita[,'2005'], 
            mode = 'lines',
            type = 'scatter3d',
            line = list(color = '#037A73',
                        width = 6)) %>%
  add_trace(x = years[1:(length(years)-1)],
            y = 23, 
            z = data3d_ita[24,], 
            mode = 'lines',
            type = 'scatter3d',
            line = list(color = '#7570b3',
                        width = 6)) %>%
  add_trace(x = years[1:(length(years)-1)],
            y = 60, 
            z = data3d_ita[61,], 
            mode = 'lines',
            type = 'scatter3d',
            line = list(color = '#7570b3',
                        width = 6)) %>%
  layout(scene = list(xaxis = list(title = 'Years', showgrid = F,
                                   linecolor = 'black',
                                   linewidth = 1.5),
                      yaxis = list(title = 'Age', showgrid = F,
                                   linecolor = 'black',
                                   linewdith = 1.5),
                      zaxis = list(title = 'log(Rate)', showgrid = F,
                                   linecolor = 'black',
                                   linewidth = 1.5),
                      camera = list(
                        center = list(x = 0, y = 0, z = -0.1),
                        eye = list(
                          x = 1.5,
                          y = 0.9,
                          z = 1.5))),
         showlegend = FALSE) %>%
  hide_colorbar() -> plot3d_ita

htmlwidgets::saveWidget(plot3d_ita,
                        file = here('output', 'ita3d.html'))

plot_ly(z=~data3d_us,
        x = years) %>%
  add_surface(alpha = 0.8,
              colors = "Spectral",
              reversescale = TRUE) %>%
  add_trace(x = rep(1970, 101), 
            y = 0:100, 
            z = data3d_us[,'1970'], 
            mode = 'lines',
            type = 'scatter3d',
            line = list(color = '#037A73',
                        width = 6)) %>%
  add_trace(x = rep(2005, 101),
            y = 0:100,
            z = data3d_us[,'2005'],
            mode = 'lines',
            type = 'scatter3d',
            line = list(color = '#037A73',
                        width = 6)) %>%
  add_trace(x = years,
            y = 23, 
            z = data3d_us[24,], 
            mode = 'lines',
            type = 'scatter3d',
            line = list(color = '#7570b3',
                        width = 6)) %>%
  add_trace(x = years,
            y = 60, 
            z = data3d_us[61,], 
            mode = 'lines',
            type = 'scatter3d',
            line = list(color = '#7570b3',
                        width = 6)) %>%
  layout(scene = list(xaxis = list(title = 'Years', showgrid = F,
                                   linecolor = 'black',
                                   linewidth = 1.5),
                      yaxis = list(title = 'Age', showgrid = F,
                                   linecolor = 'black',
                                   linewdith = 1.5),
                      zaxis = list(title = 'log(Rate)', showgrid = F,
                                   linecolor = 'black',
                                   linewidth = 1.5),
                      camera = list(
                        center = list(x = 0, y = 0, z = -0.1),
                        eye = list(
                          x = 1.5,
                          y = 0.9,
                          z = 1.5))),
         showlegend = FALSE) %>%
  hide_colorbar() -> plot3d_us

htmlwidgets::saveWidget(plot3d_us,
                        file = here('output', 'us3d.html'))



#########################################
## Data age-specific plot     ###########
#########################################

as_tibble(log((Y_ita_man + Y_ita_woman)/(N_ita_man + N_ita_woman))) %>%
  mutate(year = years[1:(length(years)-1)]) %>%
  pivot_longer(-year, names_to = 'age', values_to = 'value') %>%
  mutate_at(vars(age), as.numeric) -> data_ita

as_tibble(log((Y_us_man + Y_us_woman)/(N_us_man + N_us_woman))) %>%
  mutate(year = years) %>%
  pivot_longer(-year, names_to = 'age', values_to = 'value') %>%
  mutate_at(vars(age), as.numeric) -> data_us

data_ita %>%
  filter(year %in% c(1933,
                     1943,
                     1950,
                     1960,
                     1970,
                     1980,
                     1990,
                     2000,
                     2010,
                     2019)) %>%
  mutate(title = 'Italy (1933-2019)') %>%
  bind_rows(data_us %>%
              filter(year %in% c(1933,
                                 1943,
                                 1950,
                                 1960,
                                 1970,
                                 1980,
                                 1990,
                                 2000,
                                 2010,
                                 2020)) %>%
              mutate(title = 'United States (1933-2020)')) %>%
  ggplot(aes(x = age, y = value, color = as.factor(year))) +
  geom_line() + 
  geom_label_repel(data = . %>% 
                     filter((year == 1943) & 
                              (title == 'Italy (1933-2019)') & 
                              (age == 40) |
                              (year == 1933) &
                              (title == 'United States (1933-2020)') & 
                              (age == 40)),
                   aes(x = age, y = value, label = year),
                   nudge_y = 0.9,
                   nudge_x = 1) +
  geom_label_repel(data = . %>% 
                     filter((year == 2019) &
                              (title == 'Italy (1933-2019)') & 
                              (age == 15) |
                              (year == 2020) &
                              (title == 'United States (1933-2020)') & 
                              (age == 15)),
                   aes(x = age, y = value, label = year),
                   nudge_y = -0.5,
                   nudge_x = 15) +
  facet_grid(cols = vars(title)) +
  scale_color_brewer(palette = 'BrBG') +
  scale_y_continuous(breaks = c(-1.5,-4.5,-7.5)) +
  scale_x_continuous(n.breaks = 10) +
  guides(color = "none") +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = 'grey95'),
        plot.margin = margin(r = 0),
        axis.title = element_blank()) -> plot_ages

ggsave(plot_ages,
       filename = here('output','ita_us_ages.pdf'),
       width = 10,
       height = 4)

#########################################
## Heatmap model plot    ################
#########################################

selected <- c('U1','U6','U14')
colors_sel <- c(black = 'grey80',
                col1 = '#018571',
                col2 = '#a6611a',
                col3 = '#8c6bb1')

smoothing_tibble_adj %>%
  filter(data == 'ita_man',
         state == 'U') %>%
  left_join(age_classes, by = "weight") %>%
  mutate(label = case_when(weight == selected[1] ~ 'col1',
                           weight == selected[2] ~ 'col2',
                           weight == selected[3] ~ 'col3',
                           TRUE ~ 'black')) %>%
  mutate(title = 'Time trajectory of B-spline coefficients') %>%
  ggplot(aes(x = t, y = value, group = weight)) + 
  geom_line(aes(color = label, alpha = label), size = 0.9) +
  scale_color_manual(values = colors_sel) +
  scale_y_continuous(n.breaks = 3) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_alpha_discrete(range = c('black' = 0.3, 
                                 'col1' = 1, 
                                 'col2' = 1, 
                                 'col3' = 1)) +
  facet_grid(cols = vars(title)) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(b = 0),
        strip.background = element_rect(fill = 'grey95'),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) -> marginal_state

as_tibble(Sfull)  %>%
  mutate(age = seq(1, 100, by = 0.1)) %>%
  pivot_longer(-age, names_to = 'weight', values_to = 'value') %>%
  mutate_at(vars(value), as.numeric) %>%
  mutate(label = case_when(weight == selected[1] ~ 'col1',
                           weight == selected[2] ~ 'col2',
                           weight == selected[3] ~ 'col3',
                           TRUE ~ 'black')) %>%
  mutate(title = 'B-spline basis functions') %>%
  ggplot(aes(x = age, y = value)) +
  geom_line(aes(group = weight, color = label, size = label), alpha = 0.8) +
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = 1),
               size = 0.3,
               color = colors_spline['black']) +
  scale_size_manual(values = c("black" = 0.4, "col1" = 0.8, "col2" = 0.8, "col3" = 0.8)) +
  scale_color_manual(values = colors_sel) +
  scale_x_continuous(expand = c(0.0005,0.01), limits = c(-0.5,100)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.1)) +
  labs(x = NULL, y = NULL) +
  guides(color = "none", alpha = "none", size = "none") +
  geom_segment(aes(x = 0, y = 0, xend = 84, yend = 0), 
               color = "grey90", size = 1) +
  coord_flip() +
  facet_grid(rows = vars(title)) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = 'grey95'),
        plot.margin = margin(l = 0),
        axis.ticks.y = element_blank()) -> marginal_spline


signal_smoothing_tibble_adj %>%
  filter(data == 'ita_man',
         state == 'm') %>%
  ggplot(aes(x = t, y = age)) +
  geom_tile(aes(fill = value, color = value)) +
  scale_fill_distiller(palette = 'Spectral', direction = -1) +
  scale_color_distiller(palette = 'Spectral', direction = -1) +
  guides(color = "none") +
  labs(fill = latex2exp::TeX('$f_{t}(x)$')) +
  scale_x_continuous(expand = c(0, 0), n.breaks = 6) +
  scale_y_continuous(expand = c(0, 0), n.breaks = 8) +
  theme(axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.spacing = unit(0, "cm"),
        plot.margin = margin(r = 0, t = 0),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color = "grey90"),
        legend.title.align = 0.5,
        legend.key.size = unit(0.4, 'cm')) -> heat

(marginal_state + 
    guide_area() + 
    theme(plot.margin = margin(b = 0, r = 0))) + 
  heat + 
  marginal_spline + 
  plot_layout(ncol = 2,
              nrow  = 2, 
              guides = 'collect', 
              widths = c(9,1), 
              heights = c(1,2.5)) -> plot_model

ggsave(plot_model,
       filename = here('output','model.pdf'),
       width = 8,
       height = 5)



