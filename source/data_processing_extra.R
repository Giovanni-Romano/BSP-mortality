rm(list=ls())
require(here)
require(tidyverse)
require(magrittr)

# Data loading
fra_deaths <- read.delim(here('data','fra_deaths_period.txt'), skip=2, sep="")
fra_atrisk <- read.delim(here('data','fra_atrisk_period.txt'), skip=2, sep="")

dnk_deaths <- read.delim(here('data','dnk_deaths_period.txt'), skip=2, sep="")
dnk_atrisk <- read.delim(here('data','dnk_atrisk_period.txt'), skip=2, sep="")

cze_deaths <- read.delim(here('data','cze_deaths_period.txt'), skip=2, sep="")
cze_atrisk <- read.delim(here('data','cze_atrisk_period.txt'), skip=2, sep="")


# Years of analysis
years <- 1950:2020
ages <- 0:100

# FRA data
Y_fra_man <- fra_deaths %>%
  select(Year, Age, Male) %>%
  filter(Year %in% years,
         Age %in% ages) %>%
  pivot_wider(names_from = Age, values_from = Male) %>%
  select(-Year) %>% as.matrix() %>% round() %>%
  set_rownames(years)

Y_fra_woman <- fra_deaths %>%
  select(Year, Age, Female) %>%
  filter(Year %in% years,
         Age %in% ages) %>%
  pivot_wider(names_from = Age, values_from = Female) %>%
  select(-Year) %>% as.matrix() %>% round() %>%
  set_rownames(years)

N_fra_man <- fra_atrisk %>%
  select(Year, Age, Male) %>%
  filter(Year %in% years,
         Age %in% ages) %>%
  pivot_wider(names_from = Age, values_from = Male) %>%
  select(-Year) %>% as.matrix() %>% round() %>%
  set_rownames(years)

N_fra_woman <- fra_atrisk %>%
  select(Year, Age, Female) %>%
  filter(Year %in% years,
         Age %in% ages) %>%
  pivot_wider(names_from = Age, values_from = Female) %>%
  select(-Year) %>% as.matrix() %>% round() %>%
  set_rownames(years)

# UK data
Y_dnk_man <- dnk_deaths %>%
  select(Year, Age, Male) %>%
  filter(Year %in% years,
         Age %in% ages) %>%
  pivot_wider(names_from = Age, values_from = Male) %>%
  select(-Year) %>% as.matrix() %>% round() %>%
  set_rownames(years)
  
Y_dnk_woman <- dnk_deaths %>%
  select(Year, Age, Female) %>%
  filter(Year %in% years,
         Age %in% ages) %>%
  pivot_wider(names_from = Age, values_from = Female) %>%
  select(-Year) %>% as.matrix() %>% round() %>%
  set_rownames(years)

N_dnk_man <- dnk_atrisk %>%
  select(Year, Age, Male) %>%
  filter(Year %in% years,
         Age %in% ages) %>%
  pivot_wider(names_from = Age, values_from = Male) %>%
  select(-Year) %>% as.matrix() %>% round() %>%
  set_rownames(years)

N_dnk_woman <- dnk_atrisk %>%
  select(Year, Age, Female) %>%
  filter(Year %in% years,
         Age %in% ages) %>%
  pivot_wider(names_from = Age, values_from = Female) %>%
  select(-Year) %>% as.matrix() %>% round() %>%
  set_rownames(years)

# US data
Y_cze_man <- cze_deaths %>%
  select(Year, Age, Male) %>%
  filter(Year %in% years,
         Age %in% ages) %>%
  pivot_wider(names_from = Age, values_from = Male) %>%
  select(-Year) %>% as.matrix() %>% round() %>%
  set_rownames(years)

Y_cze_woman <- cze_deaths %>%
  select(Year, Age, Female) %>%
  filter(Year %in% years,
         Age %in% ages) %>%
  pivot_wider(names_from = Age, values_from = Female) %>%
  select(-Year) %>% as.matrix() %>% round() %>%
  set_rownames(years)

N_cze_man <- cze_atrisk %>%
  select(Year, Age, Male) %>%
  filter(Year %in% years,
         Age %in% ages) %>%
  pivot_wider(names_from = Age, values_from = Male) %>%
  select(-Year) %>% as.matrix() %>% round() %>%
  set_rownames(years)

N_cze_woman <- cze_atrisk %>%
  select(Year, Age, Female) %>%
  filter(Year %in% years,
         Age %in% ages) %>%
  pivot_wider(names_from = Age, values_from = Female) %>%
  select(-Year) %>% as.matrix() %>% round() %>%
  set_rownames(years)


# Data checks
sum(is.na(Y_fra_man))
sum(is.na(Y_fra_woman))
sum(is.na(Y_dnk_man))
sum(is.na(Y_dnk_woman))
sum(is.na(Y_cze_man))
sum(is.na(Y_cze_woman))

sum(is.na(N_fra_man))
sum(is.na(N_fra_woman))
sum(is.na(N_dnk_man))
sum(is.na(N_dnk_woman))
sum(is.na(N_cze_man))
sum(is.na(N_cze_woman))

all(dim(Y_fra_man) == dim(N_fra_man))
all(dim(Y_fra_woman) == dim(N_fra_woman))
all(dim(Y_dnk_man) == dim(N_dnk_man))
all(dim(Y_dnk_woman) == dim(N_dnk_woman))
all(dim(Y_cze_man) == dim(N_cze_man))
all(dim(Y_cze_woman) == dim(N_cze_woman))

all(Y_fra_man > 0)
all(Y_fra_woman > 0)
all(Y_dnk_man > 0)
all(Y_dnk_woman > 0)
all(Y_cze_man > 0)
all(Y_cze_woman > 0)

#### DNK and CZE have very few 0 deaths! -> changing to closest ages average
dnk_man_zero <- which(Y_dnk_man == 0, arr.ind = T)
for(it in 1:nrow(dnk_man_zero)){
  rc <- dnk_man_zero[it,]
  Y_dnk_man[rc[1],rc[2]] <- (Y_dnk_man[rc[1],rc[2]-1] + Y_dnk_man[rc[1],rc[2]+1])/2
}

dnk_woman_zero <- which(Y_dnk_woman == 0, arr.ind = T)
for(it in 1:nrow(dnk_woman_zero)){
  rc <- dnk_woman_zero[it,]
  Y_dnk_woman[rc[1],rc[2]] <- (Y_dnk_woman[rc[1],rc[2]-1] + Y_dnk_woman[rc[1],rc[2]+1])/2
}

cze_man_zero <- which(Y_cze_man == 0, arr.ind = T)
for(it in 1:nrow(cze_man_zero)){
  rc <- cze_man_zero[it,]
  Y_cze_man[rc[1],rc[2]] <- (Y_cze_man[rc[1],rc[2]-1] + 0)/2
}

cze_woman_zero <- which(Y_cze_woman == 0, arr.ind = T)
for(it in 1:nrow(cze_woman_zero)){
  rc <- cze_woman_zero[it,]
  Y_cze_woman[rc[1],rc[2]] <- (Y_cze_woman[rc[1],rc[2]-1] + Y_cze_woman[rc[1],rc[2]+1])/2
}


all(N_fra_man > 0)
all(N_fra_woman > 0)
all(N_dnk_man > 0)
all(N_dnk_woman > 0)
all(N_cze_man > 0)
all(N_cze_woman > 0)

save(Y_fra_man,
     N_fra_man,
     Y_fra_woman,
     N_fra_woman,
     Y_dnk_man,
     N_dnk_man,
     Y_dnk_woman,
     N_dnk_woman,
     Y_cze_man,
     N_cze_man,
     Y_cze_woman,
     N_cze_woman,
     years,
     ages,
     file = here('output','mortality_extra.Rdata'))
