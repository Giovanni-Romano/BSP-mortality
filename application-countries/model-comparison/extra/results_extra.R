require(here)
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
load(here('output/mortality_extra.Rdata'))
data_list <- list(fra_man = list(rates = Y_fra_man/N_fra_man),
                  fra_woman = list(rates = Y_fra_woman/N_fra_woman),
                  dnk_man = list(rates = Y_dnk_man/N_dnk_man),
                  dnk_woman = list(rates = Y_dnk_woman/N_dnk_woman),
                  cze_man = list(rates = Y_cze_man/N_cze_man),
                  cze_woman = list(rates = Y_cze_woman/N_cze_woman))
train <- 50

# Helper functions
source(here('source','helper_fun.R'))

# APC model results
load(here('output','APC_for_extra.Rdata'))

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
load(here('output','LC_for_extra.Rdata'))

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
load(here('output','CBD_for_extra.Rdata'))

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
load(here('output','PLAT_for_extra.Rdata'))

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
load(here('output','RH_for_extra.Rdata'))

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
load(here('output','HU_for_extra.Rdata'))

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
load(here('output','CP_boot_for_extra.Rdata'))

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
BSP_uq_res_age <- tibble()
for(country in c('FRA', 'DNK', 'CZE')){
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
save(list = c("BSP_uq_res_age",
              "APC_res_age",
              "LC_res_age",
              "CBD_res_age",
              "PLAT_res_age",
              "RH_res_age",
              "HU_res_age",
              "CP_res_age"),
     file = here('output','results_uq_for_extra.Rdata'))

