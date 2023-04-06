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


# # BSP model results
# BSP_res_age <- tibble()
# for(country in c('UK', 'US', 'ITA', 'SWE')){
#   load(here('output', paste(country, '_for.Rdata', sep = '')))
#   print(names(res_forward))
#   for(country_gender in names(res_forward)){
#     BSP_res_age <- bind_rows(BSP_res_age,
#                              processBSP(item = res_forward[[country_gender]],
#                                         rates = data_list[[country_gender]]$rates,
#                                         country = country) %>%
#                                mutate(data = country_gender)) %>%
#       mutate(model = 'BSP')
#   }
# }
# rm(list = 'res_forward')

## New

# # NGP model results
# NGP_res_age <- tibble()
# for(country in c('UK', 'US', 'ITA', 'SWE')){
#   load(here('output', paste(country, '_NGP_delta1_for.Rdata', sep = '')))
#   print(names(res_forward))
#   for(country_gender in names(res_forward)){
#     NGP_res_age <- bind_rows(NGP_res_age,
#                              processBSP(item = res_forward[[country_gender]],
#                                         rates = data_list[[country_gender]]$rates,
#                                         country = country) %>%
#                                mutate(data = country_gender)) %>%
#       mutate(model = 'NGP')
#   }
# }
# rm(list = 'res_forward')

# # Kalman model results
# Kalman_res_age <- tibble()
# for(country in c('UK', 'US', 'ITA', 'SWE')){
#   load(here('output', paste(country, '_Kalman_for.Rdata', sep = '')))
#   print(names(res_forward))
#   for(country_gender in names(res_forward)){
#     Kalman_res_age <- bind_rows(Kalman_res_age,
#                                 processBSP(item = res_forward[[country_gender]],
#                                            rates = data_list[[country_gender]]$rates,
#                                            country = country) %>%
#                                   mutate(data = country_gender)) %>%
#       mutate(model = 'Kalman')
#   }
# }
# rm(list = 'res_forward')

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

load(here('output','results_uq_for_extra.Rdata'))

### Result averaging over all ages, genders, and countries
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
         width95_log, inside) %>%
  group_by(model, h_ahead) %>%
  summarise(mean_abserr_log = mean(mad_log),
            median_abserr_log = median(mad_log),
            q25_abserr_log = quantile(mad_log,  probs = 0.25),
            q75_abserr_log = quantile(mad_log,  probs = 0.75),
            mean_rmse_log = mean(rmse_log),
            median_rmse_log = median(rmse_log),
            q25_rmse_log = quantile(rmse_log,  probs = 0.25),
            q75_rmse_log = quantile(rmse_log,  probs = 0.75),
            mean_width95_log = mean(width95_log),
            median_width95_log = median(width95_log),
            q25_width95_log = quantile(width95_log,  probs = 0.25),
            q75_width95_log = quantile(width95_log,  probs = 0.75),
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
  # relocate(Kalman, .before = PLAT) %>%
  # relocate(NGP, .before = PLAT) %>%
  xtable(x = ., 
         digits = 3,
         caption = 'EXTRA (FRA-DNK-CZE): Median absolute deviation of the mortality log-rate.')

print.xtable(xRes,
             file = here('output','table_median_extra.tex'),
             include.rownames = FALSE)

## RMSE
res_log <- results_tot %>%
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
  # relocate(Kalman, .before = PLAT) %>%
  # relocate(NGP, .before = PLAT) %>%
  xtable(x = ., 
         digits = 3,
         caption = 'EXTRA (FRA-DNK-CZE): Root mean squared error of the mortality log-rate.')

print.xtable(xRes,
             file = here('output','table_rmse_extra.tex'),
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
  # relocate(Kalman, .before = PLAT) %>%
  # relocate(NGP, .before = PLAT) %>%
  xtable(x = ., 
         digits = 3,
         caption = 'EXTRA (FRA-DNK-CZE): Coverage of 95\\% prediction interval')

print.xtable(xRes,
             file = here('output','table_uq_extra.tex'),
             include.rownames = FALSE)

## Width 95
res_log <- results_tot %>%
  select(mean_width95_log, model, h_ahead) %>%
  pivot_wider(names_from = "model", values_from = 'mean_width95_log')

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
  # relocate(Kalman, .before = PLAT) %>%
  # relocate(NGP, .before = PLAT) %>%
  xtable(x = ., 
         digits = 3,
         caption = 'EXTRA (FRA-DNK-CZE): Average width of 95\\% prediction interval')

print.xtable(xRes,
             file = here('output','table_width_extra.tex'),
             include.rownames = FALSE)



