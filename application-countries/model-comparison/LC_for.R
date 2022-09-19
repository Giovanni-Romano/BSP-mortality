require(here)
require(tidyverse)
require(StMoMo)
require(doParallel)

rm(list=ls())

set.seed(3015)

load(here('output','mortality.Rdata'))
output_collector <- list()
output_collector <- append(output_collector, list(info=Sys.time()))

# LC model via StMoMo package
LC <- lc(link = 'log') 

# Set up minimum training window, how many step-ahead to predict
train <- 58
h_step <- 10
output_collector <- append(output_collector, list(train = train, h_step = h_step))

# Function to fit  and forecast
fit_sim <- function(t, train, mod, h, stmomo_data, years, gc.order = NULL, nsim = 10){
  wxt_osa <- genWeightMat(ages = stmomo_data$ages,
                          years = years[1:t], clip = 3)
  
  err2 <- err <- try(log('a'), silent = TRUE) # random error
  while(inherits(err, "try-error") | inherits(err2, "try-error")){ # keep repeating if there is an error
    err <- try(fit <- fit(mod,
                          data = stmomo_data,
                          wxt = wxt_osa,
                          gc.order = gc.order,
                          years = years,
                          years.fit = years[1:t],
                          verbose = FALSE))
    print('fit done')
    err2 <- try(sim <- simulate(fit, h = h, nsim = nsim))
  }
  return(sim)
}

# Function to set up the rolling window: minimum window given by train, then 
# increasing the training set by 1 each time
rolling <- function(cg){
  country <- sub("_.*", "", cg)
  gender <- sub(".*_", "", cg)
  print(paste('Doing', country, gender))
  Y <- eval(parse(text = paste('Y', country, gender, sep = '_')))
  N <- eval(parse(text = paste('N', country, gender, sep = '_')))
  
  years <- years[1:nrow(Y)]
  
  stmomo_data <- structure(list(Dxt = t(Y), 
                                Ext = t(N), 
                                ages = ages, 
                                years = years, 
                                type = 'central', 
                                series = ifelse(gender == 'man', 'male', 'female'),  # whether "male", "female" or "total"
                                label = country), 
                           class = "StMoMoData")
  
  wxt <- genWeightMat(ages = stmomo_data$ages, 
                      years = stmomo_data$years, 
                      clip = 3)
  
  registerDoParallel(cores = 50)
  cl <- makeCluster(50, type = "FORK")
  LCsim_ <- parLapply(cl,
                      train:length(years),
                      . %>%
                        fit_sim(.,
                                train = train,
                                mod = LC,
                                h = h_step,
                                stmomo_data = stmomo_data,
                                years = years,
                                nsim = 100))
  stopCluster(cl)
  
  LCsim <- LCsim_ %>%
    modify(. %>% pluck('rates'))
  return(LCsim)
}

# Set up parameter with country-gender combination
cg <- expand.grid(country = c('uk','us','swe','ita'),
                  gender = c('man', 'woman')) %>%
  mutate(cg = paste(country, gender, sep = '_')) %>%
  pull(cg)

# Fit and forecast via lapply
res_forward <- lapply(cg, 
                      rolling)
names(res_forward) <- cg

output_collector <- append(output_collector, list(warnings = warnings()))

save(list = c('output_collector',
              'res_forward'),
     file = here('output','LC_for.Rdata'))
