require(tidyverse)
require(ggthemes)
theme_set(theme_light() + theme(panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                strip.background = element_rect(color='gray',fill='white'),
                                strip.text.x = element_text(color='#525252'),
                                strip.text.y = element_text(color='#525252')))
require(xtable)

rm(list = ls())

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
load(here('output','CP_for.Rdata'))

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


# BSP model results
BSP_res_age <- tibble()
for(country in c('UK', 'US', 'ITA', 'SWE')){
  load(here('output', paste(country, '_for.Rdata', sep = '')))
  print(names(res_forward))
  for(country_gender in names(res_forward)){
    BSP_res_age <- bind_rows(BSP_res_age,
                             processBSP(item = res_forward[[country_gender]],
                                        rates = data_list[[country_gender]]$rates,
                                        country = country) %>%
                               mutate(data = country_gender)) %>%
      mutate(model = 'BSP')
  }
}
rm(list = 'res_forward')



# Results
save(list = c("BSP_res_age",
              "APC_res_age",
              "LC_res_age",
              "CBD_res_age",
              "PLAT_res_age",
              "RH_res_age",
              "HU_res_age",
              "CP_res_age"),
     file = here('output','results_for.Rdata'))

# load(here('results_for.Rdata'))

#### Result averaging over all ages, genders, and countries

bind_rows(BSP_res_age,
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
         rmse_log, mad, mad_log) %>%
  group_by(model, h_ahead) %>%
  summarise(mean_abserr_log = mean(mad_log),
            median_abserr_log = median(mad_log)) %>%
  ungroup() -> results_tot

#### Overall on log scale ####
res_log <- results_tot %>%
  select(median_abserr_log, model, h_ahead) %>%
  pivot_wider(names_from = "model", values_from = 'median_abserr_log')

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
         caption = 'Median absolute deviation of the mortality log-rate.')

print.xtable(xRes,
             file = here('output','table.tex'),
             include.rownames = FALSE)

##
res_log <- results_tot %>%
  select(mean_abserr_log, model, h_ahead) %>%
  pivot_wider(names_from = "model", values_from = 'mean_abserr_log')

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
         caption = 'Mean absolute deviation of the mortality log-rate.')

print.xtable(xRes,
             file = here('output','table.tex'),
             append = TRUE,
             include.rownames = FALSE)







