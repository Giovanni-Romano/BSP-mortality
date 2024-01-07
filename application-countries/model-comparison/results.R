require(here)
require(tidyverse)
require(ggthemes)
theme_set(theme_light() + theme(panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                strip.background = element_rect(color='gray',fill='white'),
                                strip.text.x = element_text(color='#525252'),
                                strip.text.y = element_text(color='#525252')))
require(xtable)
require(patchwork)

rm(list = ls())

dem.env <- getNamespace("demography") # requires package "demography"

set.seed(5674)

# Data
load(here('output/mortality.Rdata'))
data_list <- list(uk_man = list(rates = Y_uk_man/N_uk_man),
                  uk_woman = list(rates = Y_uk_woman/N_uk_woman),
                  us_man = list(rates = Y_us_man/N_us_man),
                  us_woman = list(rates = Y_us_woman/N_us_woman),
                  swe_man = list(rates = Y_swe_man/N_swe_man),
                  swe_woman = list(rates = Y_swe_woman/N_swe_woman),
                  ita_man = list(rates = Y_ita_man/N_ita_man),
                  ita_woman = list(rates = Y_ita_woman/N_ita_woman))
years <- 1933:2020
ages <- 0:100
train <- 58

# Helper functions
source(here('source','helper_fun.R'))

# APC model results
load(here('output','APC_for.Rdata'))

APC_res_age <- tibble()
for(country_gender in names(res_forward)){
  print(paste('APC', country_gender))
  country <- toupper(sub("_.*", "", country_gender))
  APC_res_age <- bind_rows(APC_res_age,
                           processStMoMo(item = res_forward[[country_gender]], 
                                         model = 'APC', 
                                         rates = data_list[[country_gender]]$rates,
                                         country = country) %>%
                             mutate(data = country_gender))
}
rm(list = 'res_forward')

# LC model results
load(here('output','LC_for.Rdata'))

LC_res_age <- tibble()
for(country_gender in names(res_forward)){
  print(paste('LC', country_gender))
  country <- toupper(sub("_.*", "", country_gender))
  LC_res_age <- bind_rows(LC_res_age,
                           processStMoMo(item = res_forward[[country_gender]], 
                                         model = 'LC', 
                                         rates = data_list[[country_gender]]$rates,
                                         country = country) %>%
                             mutate(data = country_gender))
}

rm(list = 'res_forward')

# CBD model results
load(here('output','CBD_for.Rdata'))

CBD_res_age <- tibble()
for(country_gender in names(res_forward)){
  print(paste('CBD', country_gender))
  country <- toupper(sub("_.*", "", country_gender))
  CBD_res_age <- bind_rows(CBD_res_age,
                          processStMoMo(item = res_forward[[country_gender]], 
                                        model = 'CBD', 
                                        rates = data_list[[country_gender]]$rates,
                                        country = country) %>%
                            mutate(data = country_gender))
}
rm(list = 'res_forward')

# PLAT model results
load(here('output','PLAT_for.Rdata'))

PLAT_res_age <- tibble()
for(country_gender in names(res_forward)){
  print(paste('PLAT', country_gender))
  country <- toupper(sub("_.*", "", country_gender))
  PLAT_res_age <- bind_rows(PLAT_res_age,
                          processStMoMo(item = res_forward[[country_gender]], 
                                        model = 'PLAT', 
                                        rates = data_list[[country_gender]]$rates,
                                        country = country) %>%
                            mutate(data = country_gender))
}
rm(list = 'res_forward')

# RH model results
load(here('output','RH_for.Rdata'))

RH_res_age <- tibble()
for(country_gender in names(res_forward)){
  print(paste('RH', country_gender))
  country <- toupper(sub("_.*", "", country_gender))
  RH_res_age <- bind_rows(RH_res_age,
                            processStMoMo(item = res_forward[[country_gender]], 
                                          model = 'RH', 
                                          rates = data_list[[country_gender]]$rates,
                                          country = country) %>%
                              mutate(data = country_gender))
}
rm(list = 'res_forward')

# HU model results
load(here('output','HU_for.Rdata'))

HU_res_age <- tibble()
for(country_gender in names(res_forward)){
  print(paste('HU', country_gender))
  country <- toupper(sub("_.*", "", country_gender))
  HU_res_age <- bind_rows(HU_res_age,
                      processHU(item = res_forward[[country_gender]],
                                rates = data_list[[country_gender]]$rates,
                                country = country) %>%
                        mutate(data = country_gender))
}
rm(list = 'res_forward')

# CP model results
load(here('output','CP_boot_for.Rdata'))

CP_res_age <- tibble()
for(country_gender in names(res_forward)){
  print(paste('CP', country_gender))
  country <- toupper(sub("_.*", "", country_gender))
  CP_res_age <- bind_rows(CP_res_age,
                          processCP(item = res_forward[[country_gender]],
                                    rates = data_list[[country_gender]]$rates,
                                    country = country) %>%
                            mutate(data = country_gender))
}
rm(list = 'res_forward')

# NGP model results
NGP_res_age <- tibble()
for(country in c('UK', 'US', 'ITA', 'SWE')){
  load(here('output', paste(country, '_NGP_delta1_for.Rdata', sep = '')))
  print(names(res_forward))
  for(country_gender in names(res_forward)){
    NGP_res_age <- bind_rows(NGP_res_age,
                             processBSP(item = res_forward[[country_gender]],
                                        rates = data_list[[country_gender]]$rates,
                                        country = country) %>%
                               mutate(data = country_gender)) %>%
      mutate(model = 'NGP')
  }
}
rm(list = 'res_forward')

# Kalman model results
Kalman_res_age <- tibble()
for(country in c('UK', 'US', 'ITA', 'SWE')){
  load(here('output', paste(country, '_Kalman_for.Rdata', sep = '')))
  print(names(res_forward))
  for(country_gender in names(res_forward)){
    Kalman_res_age <- bind_rows(Kalman_res_age,
                                processBSP(item = res_forward[[country_gender]],
                                           rates = data_list[[country_gender]]$rates,
                                           country = country) %>%
                                  mutate(data = country_gender)) %>%
      mutate(model = 'Kalman')
  }
}
rm(list = 'res_forward')

# BSP model results
BSP_uq_res_age <- tibble()
for(country in c('UK', 'US', 'ITA', 'SWE')){
  load(here('output', paste(country, '_uq_for.Rdata', sep = '')))
  print(names(res_forward))
  for(country_gender in names(res_forward)){
    BSP_uq_res_age <- bind_rows(BSP_uq_res_age,
                                processBSP(item = res_forward[[country_gender]],
                                           rates = data_list[[country_gender]]$rates,
                                           country = country) %>%
                                  mutate(data = country_gender)) %>%
      mutate(model = 'BSP')
  }
}
rm(list = 'res_forward')

# Results
# save(list = c("BSP_uq_res_age",
#               "APC_res_age",
#               "LC_res_age",
#               "CBD_res_age",
#               "PLAT_res_age",
#               "RH_res_age",
#               "NGP_res_age",
#               "Kalman_res_age",
#               "HU_res_age",
#               "CP_res_age"),
#      file = here('output','results_uq_for.Rdata'))

load(here('output','results_uq_for.Rdata'))

### Result averaging over all ages, genders, and countries
bind_rows(BSP_uq_res_age,
          NGP_res_age,
          Kalman_res_age,
          APC_res_age,
          LC_res_age,
          CBD_res_age,
          PLAT_res_age,
          RH_res_age,
          HU_res_age,
          CP_res_age) %>%
  filter(time_fit >= 1990,
         time_fit <= 2010,
         time_pred <= 2020,
         (data == 'uk_man') | (data == 'uk_woman') |
           (data == 'us_man') | (data == 'us_woman') |
           (data == 'ita_man' & time_fit <= 2009 & time_pred <= 2019) | 
           (data == 'ita_woman' & time_fit <= 2009 & time_pred <= 2019) |
           (data == 'swe_man') | (data == 'swe_woman')) %>%
  select(time_fit, age, h_ahead, model,
         data, time_pred, rmse,
         rmse_log, mad, mad_log,
         width95_log, inside) %>%
  group_by(model, h_ahead) %>%
  summarise(median_abserr_log = median(mad_log),
            q25_abserr_log = quantile(mad_log,  probs = 0.25),
            q75_abserr_log = quantile(mad_log,  probs = 0.75),
            coverage = mean(inside)) %>%
  ungroup() -> results_tot


#### Overall on log scale ####
res_log <- results_tot %>%
  select(median_abserr_log, model, h_ahead) %>%
  pivot_wider(names_from = "model", values_from = 'median_abserr_log')

res_quantiles_log <- results_tot %>%
  filter(model %in% c("BSP", "CP", "HU")) %>%
  select(q25_abserr_log, model, h_ahead) %>%
  pivot_wider(names_from = "model", values_from = 'q25_abserr_log') %>%
  left_join(results_tot %>%
              filter(model %in% c("BSP", "CP", "HU")) %>%
              select(q75_abserr_log, model, h_ahead) %>%
              pivot_wider(names_from = "model", values_from = 'q75_abserr_log'),
            by = c("h_ahead" = "h_ahead"),
            suffix = c(".25",".75"))

options(xtable.floating = TRUE)
options(xtable.timestamp = "")

xRes <- res_log %>%
  left_join(res_quantiles_log,
            by = c("h_ahead" = "h_ahead")) %>%
  mutate_at(vars(BSP), ~ paste(round(., digits = 3), " [", 
                               round(BSP.25, digits = 2), ",", 
                               round(BSP.75, digits = 2), 
                               "]", sep = "")) %>%
  mutate_at(vars(CP), ~ paste(round(., digits = 3), " [", 
                               round(CP.25, digits = 2), ",", 
                               round(CP.75, digits = 2), 
                               "]", sep = "")) %>%
  mutate_at(vars(HU), ~ paste(round(., digits = 3), " [", 
                               round(HU.25, digits = 2), ",", 
                               round(HU.75, digits = 2), 
                               "]", sep = "")) %>%
  select(-c("BSP.25", "BSP.75",
            "CP.25", "CP.75",
            "HU.25", "HU.75")) %>%
  mutate_at(vars(h_ahead), as.integer) %>%
  rename(`Steps ahead` = h_ahead) %>%
  relocate(BSP, .before = APC) %>%
  relocate(CP, .before = APC) %>%
  relocate(HU, .before = APC) %>%
  relocate(PLAT, .before = APC) %>%
  relocate(RH, .before = APC) %>%
  relocate(LC, .before = CBD) %>%
  relocate(Kalman, .before = PLAT) %>%
  relocate(NGP, .before = PLAT) %>%
  xtable(x = ., 
         digits = 3,
         caption = 'Median absolute deviation of the mortality log-rate.')

print.xtable(xRes,
             file = here('output','table_median.tex'),
             include.rownames = FALSE)


## Coverage 95
res_log <- results_tot %>%
  select(coverage, model, h_ahead) %>%
  pivot_wider(names_from = "model", values_from = 'coverage')

options(xtable.floating = TRUE)
options(xtable.timestamp = "")

xRes <- res_log %>%
  mutate_at(vars(h_ahead), as.integer) %>%
  rename(`Steps ahead` = h_ahead) %>%
  relocate(BSP, .before = APC) %>%
  relocate(CP, .before = APC) %>%
  relocate(HU, .before = APC) %>%
  relocate(PLAT, .before = APC) %>%
  relocate(RH, .before = APC) %>%
  relocate(LC, .before = CBD) %>%
  relocate(Kalman, .before = PLAT) %>%
  relocate(NGP, .before = PLAT) %>%
  xtable(x = ., 
         digits = 3,
         caption = 'Coverage of 95\\% prediction interval')

print.xtable(xRes,
             file = here('output','table_uq.tex'),
             include.rownames = FALSE)


# Life expectancy
compute_lifexp <- function(mx){
  dem.env$lt(mx,0,agegroup = 1,sex = 'woman')$ex # Sex does not matter
}

bind_rows(BSP_uq_res_age,
          NGP_res_age,
          Kalman_res_age,
          APC_res_age,
          LC_res_age,
          CBD_res_age,
          PLAT_res_age,
          RH_res_age,
          HU_res_age,
          CP_res_age) %>%
  filter(time_fit >= 1990,
         time_fit <= 2010,
         time_pred <= 2020,
         (data == 'uk_man') | (data == 'uk_woman') |
           (data == 'us_man') | (data == 'us_woman') |
           (data == 'ita_man' & time_fit <= 2009 & time_pred <= 2019) |
           (data == 'ita_woman' & time_fit <= 2009 & time_pred <= 2019) |
           (data == 'swe_man') | (data == 'swe_woman')) %>%
  select(time_fit, age, h_ahead, model,
         data, time_pred, fit, rate.obs) %>%
  group_by(data, model, h_ahead, time_fit) %>%
  arrange(age) %>%
  mutate(exp.fit = compute_lifexp(fit),
         exp.obs = compute_lifexp(rate.obs)) %>%
  ungroup() %>%
  filter(age == 0) %>%
  group_by(model, h_ahead) %>%
  summarise(exp.median_abserr = median(abs(exp.fit - exp.obs)),
            exp.q75_abserr = quantile(abs(exp.fit - exp.obs), probs = 0.75),
            exp.q25_abserr = quantile(abs(exp.fit - exp.obs), probs = 0.25)) %>%
  ungroup() -> res_exp_age0


# MAD for life exp at age 0
res_mad_exp0 <- res_exp_age0 %>%
  select(exp.median_abserr, model, h_ahead) %>%
  pivot_wider(names_from = "model", values_from = 'exp.median_abserr')

res_quantiles_exp0 <- res_exp_age0 %>%
  filter(model %in% c("BSP", "CP", "HU")) %>%
  select(exp.q25_abserr, model, h_ahead) %>%
  pivot_wider(names_from = "model", values_from = 'exp.q25_abserr') %>%
  left_join(res_exp_age0 %>%
              filter(model %in% c("BSP", "CP", "HU")) %>%
              select(exp.q75_abserr, model, h_ahead) %>%
              pivot_wider(names_from = "model", values_from = 'exp.q75_abserr'),
            by = c("h_ahead" = "h_ahead"),
            suffix = c(".25",".75"))

options(xtable.floating = TRUE)
options(xtable.timestamp = "")

xRes <- res_mad_exp0 %>%
  left_join(res_quantiles_exp0,
            by = c("h_ahead" = "h_ahead")) %>%
  mutate_at(vars(BSP), ~ paste(round(., digits = 2), " [",
                               round(BSP.25, digits = 2), ",",
                               round(BSP.75, digits = 2),
                               "]", sep = "")) %>%
  mutate_at(vars(CP), ~ paste(round(., digits = 2), " [",
                              round(CP.25, digits = 2), ",",
                              round(CP.75, digits = 2),
                              "]", sep = "")) %>%
  mutate_at(vars(HU), ~ paste(round(., digits = 2), " [",
                              round(HU.25, digits = 2), ",",
                              round(HU.75, digits = 2),
                              "]", sep = "")) %>%
  select(-c("BSP.25", "BSP.75",
            "CP.25", "CP.75",
            "HU.25", "HU.75")) %>%
  mutate_at(vars(h_ahead), as.integer) %>%
  rename(`Steps ahead` = h_ahead) %>%
  relocate(BSP, .before = APC) %>%
  relocate(CP, .before = APC) %>%
  relocate(HU, .before = APC) %>%
  relocate(PLAT, .before = APC) %>%
  relocate(RH, .before = APC) %>%
  relocate(LC, .before = CBD) %>%
  relocate(Kalman, .before = PLAT) %>%
  relocate(NGP, .before = PLAT) %>%
  select(-c(Kalman,NGP)) %>%
  xtable(x = .,
         digits = 2,
         caption = 'Median absolute deviation of the life expectancy at birth.')

print.xtable(xRes,
             file = here('output','table_median_exp0.tex'),
             include.rownames = FALSE)


####### ALL
### ITA US UK SWE FRA DNK CZE
bind_rows(BSP_uq_res_age,
          APC_res_age,
          LC_res_age,
          CBD_res_age,
          PLAT_res_age,
          RH_res_age,
          HU_res_age,
          CP_res_age) %>%
  filter(time_fit >= 1990,
         time_fit <= 2010,
         time_pred <= 2020,
         (data == 'uk_man') | (data == 'uk_woman') |
           (data == 'us_man') | (data == 'us_woman') |
           (data == 'ita_man' & time_fit <= 2009 & time_pred <= 2019) | 
           (data == 'ita_woman' & time_fit <= 2009 & time_pred <= 2019) |
           (data == 'swe_man') | (data == 'swe_woman')) %>%
  select(time_fit, age, h_ahead, model,
         data, time_pred, rmse,
         rmse_log, mad, mad_log,
         width95_log, inside) -> res_base

load(here('output','results_uq_for_extra.Rdata'))
bind_rows(BSP_uq_res_age,
          APC_res_age,
          LC_res_age,
          CBD_res_age,
          PLAT_res_age,
          RH_res_age,
          CP_res_age,
          HU_res_age) %>%
  filter(time_pred <= 2020) %>%
  select(time_fit, age, h_ahead, model,
         data, time_pred, rmse,
         rmse_log, mad, mad_log,
         width95_log, inside) -> res_extra

res_base %>%
  bind_rows(res_extra) %>%
  group_by(model, h_ahead) %>%
  summarise(median_abserr_log = median(mad_log),
            q25_abserr_log = quantile(mad_log,  probs = 0.25),
            q75_abserr_log = quantile(mad_log,  probs = 0.75),
            mean_rmse_log = mean(rmse_log),
            q25_rmse_log = quantile(rmse_log,  probs = 0.25),
            q75_rmse_log = quantile(rmse_log,  probs = 0.75)) %>%
  ungroup() -> results_all

#### Overall on log scale ####
res_log <- results_all %>%
  select(median_abserr_log, model, h_ahead) %>%
  pivot_wider(names_from = "model", values_from = 'median_abserr_log')

res_quantiles_log <- results_all %>%
  filter(model %in% c("BSP", "CP", "HU")) %>%
  select(q25_abserr_log, model, h_ahead) %>%
  pivot_wider(names_from = "model", values_from = 'q25_abserr_log') %>%
  left_join(results_all %>%
              filter(model %in% c("BSP", "CP", "HU")) %>%
              select(q75_abserr_log, model, h_ahead) %>%
              pivot_wider(names_from = "model", values_from = 'q75_abserr_log'),
            by = c("h_ahead" = "h_ahead"),
            suffix = c(".25",".75"))

options(xtable.floating = TRUE)
options(xtable.timestamp = "")

xRes <- res_log %>%
  left_join(res_quantiles_log,
            by = c("h_ahead" = "h_ahead")) %>%
  mutate_at(vars(BSP), ~ paste(round(., digits = 3), " [", 
                               round(BSP.25, digits = 2), ",", 
                               round(BSP.75, digits = 2), 
                               "]", sep = "")) %>%
  mutate_at(vars(CP), ~ paste(round(., digits = 3), " [", 
                              round(CP.25, digits = 2), ",", 
                              round(CP.75, digits = 2), 
                              "]", sep = "")) %>%
  mutate_at(vars(HU), ~ paste(round(., digits = 3), " [", 
                              round(HU.25, digits = 2), ",", 
                              round(HU.75, digits = 2), 
                              "]", sep = "")) %>%
  select(-c("BSP.25", "BSP.75",
            "CP.25", "CP.75",
            "HU.25", "HU.75")) %>%
  mutate_at(vars(h_ahead), as.integer) %>%
  rename(`Steps ahead` = h_ahead) %>%
  relocate(BSP, .before = APC) %>%
  relocate(CP, .before = APC) %>%
  relocate(HU, .before = APC) %>%
  relocate(PLAT, .before = APC) %>%
  relocate(RH, .before = APC) %>%
  relocate(LC, .before = CBD) %>%
  xtable(x = ., 
         digits = 3,
         caption = 'ALL: Median absolute deviation of the mortality log-rate.')

print.xtable(xRes,
             file = here('output','table_median_all.tex'),
             include.rownames = FALSE)

## RMSE
res_log <- results_all %>%
  select(mean_rmse_log, model, h_ahead) %>%
  pivot_wider(names_from = "model", values_from = 'mean_rmse_log')

options(xtable.floating = TRUE)
options(xtable.timestamp = "")

xRes <- res_log %>%
  mutate_at(vars(h_ahead), as.integer) %>%
  rename(`Steps ahead` = h_ahead) %>%
  relocate(BSP, .before = APC) %>%
  relocate(CP, .before = APC) %>%
  relocate(HU, .before = APC) %>%
  relocate(PLAT, .before = APC) %>%
  relocate(RH, .before = APC) %>%
  relocate(LC, .before = CBD) %>%
  xtable(x = ., 
         digits = 3,
         caption = 'ALL: RMSE of the mortality log-rate.')

print.xtable(xRes,
             file = here('output','table_rmse_all.tex'),
             include.rownames = FALSE)
