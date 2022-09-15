rm(list=ls())
require(tidyverse)
require(magrittr)

# Data loading
ita_deaths <- read.delim(here('data','ita_deaths_period.txt'), skip=2, sep="")
ita_atrisk <- read.delim(here('data','ita_atrisk_period.txt'), skip=2, sep="")

uk_deaths <- read.delim(here('data','uk_deaths_period.txt'), skip=2, sep="")
uk_atrisk <- read.delim(here('data','uk_atrisk_period.txt'), skip=2, sep="")

us_deaths <- read.delim(here('data','us_deaths_period.txt'), skip=2, sep="")
us_atrisk <- read.delim(here('data','us_atrisk_period.txt'), skip=2, sep="")

swe_deaths <- read.delim(here('data','swe_deaths_period.txt'), skip=2, sep="")
swe_atrisk <- read.delim(here('data','swe_atrisk_period.txt'), skip=2, sep="")


# Years of analysis
years <- 1933:2020
ages <- 0:100

# ITA data
Y_ita_man <- ita_deaths %>%
  select(Year, Age, Male) %>%
  filter(Year %in% years,
         Age %in% ages) %>%
  pivot_wider(names_from = Age, values_from = Male) %>%
  select(-Year) %>% as.matrix() %>% round() %>%
  set_rownames(years[1:(length(years)-1)])

Y_ita_woman <- ita_deaths %>%
  select(Year, Age, Female) %>%
  filter(Year %in% years,
         Age %in% ages) %>%
  pivot_wider(names_from = Age, values_from = Female) %>%
  select(-Year) %>% as.matrix() %>% round() %>%
  set_rownames(years[1:(length(years)-1)])

N_ita_man <- ita_atrisk %>%
  select(Year, Age, Male) %>%
  filter(Year %in% years,
         Age %in% ages) %>%
  pivot_wider(names_from = Age, values_from = Male) %>%
  select(-Year) %>% as.matrix() %>% round() %>%
  set_rownames(years[1:(length(years)-1)])

N_ita_woman <- ita_atrisk %>%
  select(Year, Age, Female) %>%
  filter(Year %in% years,
         Age %in% ages) %>%
  pivot_wider(names_from = Age, values_from = Female) %>%
  select(-Year) %>% as.matrix() %>% round() %>%
  set_rownames(years[1:(length(years)-1)])

# UK data
Y_uk_man <- uk_deaths %>%
  select(Year, Age, Male) %>%
  filter(Year %in% years,
         Age %in% ages) %>%
  pivot_wider(names_from = Age, values_from = Male) %>%
  select(-Year) %>% as.matrix() %>% round() %>%
  set_rownames(years)
  
Y_uk_woman <- uk_deaths %>%
  select(Year, Age, Female) %>%
  filter(Year %in% years,
         Age %in% ages) %>%
  pivot_wider(names_from = Age, values_from = Female) %>%
  select(-Year) %>% as.matrix() %>% round() %>%
  set_rownames(years)

N_uk_man <- uk_atrisk %>%
  select(Year, Age, Male) %>%
  filter(Year %in% years,
         Age %in% ages) %>%
  pivot_wider(names_from = Age, values_from = Male) %>%
  select(-Year) %>% as.matrix() %>% round() %>%
  set_rownames(years)

N_uk_woman <- uk_atrisk %>%
  select(Year, Age, Female) %>%
  filter(Year %in% years,
         Age %in% ages) %>%
  pivot_wider(names_from = Age, values_from = Female) %>%
  select(-Year) %>% as.matrix() %>% round() %>%
  set_rownames(years)

# US data
Y_us_man <- us_deaths %>%
  select(Year, Age, Male) %>%
  filter(Year %in% years,
         Age %in% ages) %>%
  pivot_wider(names_from = Age, values_from = Male) %>%
  select(-Year) %>% as.matrix() %>% round() %>%
  set_rownames(years)

Y_us_woman <- us_deaths %>%
  select(Year, Age, Female) %>%
  filter(Year %in% years,
         Age %in% ages) %>%
  pivot_wider(names_from = Age, values_from = Female) %>%
  select(-Year) %>% as.matrix() %>% round() %>%
  set_rownames(years)

N_us_man <- us_atrisk %>%
  select(Year, Age, Male) %>%
  filter(Year %in% years,
         Age %in% ages) %>%
  pivot_wider(names_from = Age, values_from = Male) %>%
  select(-Year) %>% as.matrix() %>% round() %>%
  set_rownames(years)

N_us_woman <- us_atrisk %>%
  select(Year, Age, Female) %>%
  filter(Year %in% years,
         Age %in% ages) %>%
  pivot_wider(names_from = Age, values_from = Female) %>%
  select(-Year) %>% as.matrix() %>% round() %>%
  set_rownames(years)

# SWE data
Y_swe_man <- swe_deaths %>%
  select(Year, Age, Male) %>%
  filter(Year %in% years,
         Age %in% ages) %>%
  pivot_wider(names_from = Age, values_from = Male) %>%
  select(-Year) %>% as.matrix() %>% round() %>%
  set_rownames(years)

Y_swe_woman <- swe_deaths %>%
  select(Year, Age, Female) %>%
  filter(Year %in% years,
         Age %in% ages) %>%
  pivot_wider(names_from = Age, values_from = Female) %>%
  select(-Year) %>% as.matrix() %>% round() %>%
  set_rownames(years)

N_swe_man <- swe_atrisk %>%
  select(Year, Age, Male) %>%
  filter(Year %in% years,
         Age %in% ages) %>%
  pivot_wider(names_from = Age, values_from = Male) %>%
  select(-Year) %>% as.matrix() %>% round() %>%
  set_rownames(years)

N_swe_woman <- swe_atrisk %>%
  select(Year, Age, Female) %>%
  filter(Year %in% years,
         Age %in% ages) %>%
  pivot_wider(names_from = Age, values_from = Female) %>%
  select(-Year) %>% as.matrix() %>% round() %>%
  set_rownames(years)

# Data checks
sum(is.na(Y_ita_man))
sum(is.na(Y_ita_woman))
sum(is.na(Y_uk_man))
sum(is.na(Y_uk_woman))
sum(is.na(Y_us_man))
sum(is.na(Y_us_woman))
sum(is.na(Y_swe_man))
sum(is.na(Y_swe_woman))

sum(is.na(N_ita_man))
sum(is.na(N_ita_woman))
sum(is.na(N_uk_man))
sum(is.na(N_uk_woman))
sum(is.na(N_us_man))
sum(is.na(N_us_woman))
sum(is.na(N_swe_man))
sum(is.na(N_swe_woman))

all(dim(Y_ita_man) == dim(N_ita_man))
all(dim(Y_ita_woman) == dim(N_ita_woman))
all(dim(Y_uk_man) == dim(N_uk_man))
all(dim(Y_uk_woman) == dim(N_uk_woman))
all(dim(Y_us_man) == dim(N_us_man))
all(dim(Y_us_woman)==dim(N_us_woman))
all(dim(Y_swe_man) == dim(N_swe_man))
all(dim(Y_swe_woman) == dim(N_swe_woman))

all(Y_ita_man > 0)
all(Y_ita_woman > 0)
all(Y_uk_man > 0)
all(Y_uk_woman > 0)
all(Y_us_man > 0)
all(Y_us_woman > 0)
all(Y_swe_man > 0)
all(Y_swe_woman > 0)

#### Sweden has very few 0 deaths! -> changing to closest ages average
swe_man_zero <- which(Y_swe_man == 0, arr.ind = T)
for(it in 1:nrow(swe_man_zero)){
  rc <- swe_man_zero[it,]
  Y_swe_man[rc[1],rc[2]] <- (Y_swe_man[rc[1],rc[2]-1] + Y_swe_man[rc[1],rc[2]+1])/2
}

swe_woman_zero <- which(Y_swe_woman == 0, arr.ind = T)
for(it in 1:nrow(swe_woman_zero)){
  rc <- swe_woman_zero[it,]
  Y_swe_woman[rc[1],rc[2]] <- (Y_swe_woman[rc[1],rc[2]-1] + Y_swe_woman[rc[1],rc[2]+1])/2
}


all(N_ita_man > 0)
all(N_ita_woman > 0)
all(N_uk_man > 0)
all(N_uk_woman > 0)
all(N_us_man > 0)
all(N_us_woman > 0)
all(N_swe_man > 0)
all(N_swe_woman > 0)

save(Y_ita_man,
     N_ita_man,
     Y_ita_woman,
     N_ita_woman,
     Y_uk_man,
     N_uk_man,
     Y_uk_woman,
     N_uk_woman,
     Y_us_man,
     N_us_man,
     Y_us_woman,
     N_us_woman,
     Y_swe_man,
     N_swe_man,
     Y_swe_woman,
     N_swe_woman,
     years,
     ages,
     file = here('output','mortality.Rdata'))
