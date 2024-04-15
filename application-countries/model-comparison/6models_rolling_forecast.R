setwd(this.path::here())
require(here)
require(tidyverse)
require(KFAS)
require(parallel)

rm(list=ls())

set.seed(5424)
RNGkind("L'Ecuyer-CMRG")

output_collector <- list()
output_collector <- append(output_collector, list(info=Sys.time()))

## Parameters
# n_for: last *n_for* years to use to build forecasting model
# train: minimum training size
# h_step: how many steps-ahead to forecast
n_for <- 25  
train <- 58  
h_step <- 10 
output_collector <- append(output_collector, list(train = train, 
                                                  h_step = h_step))

# Mortality data
load(here('output','mortality.Rdata'))

years <- years[1:nrow(Y_ita_man)]

# Model definition
source(here('source','BSP.R'))
source(here('source','setup.R'))

output_collector <- append(output_collector, list(delta = delta,
                                                  age_knots = age_knots,
                                                  rep = rep,
                                                  lik = 'gaussian',
                                                  kernel = matern_kernel,
                                                  n_for = n_for))

# Loading forecast and forward-rolling functions
source(here('source','BSP_forecast.R'))

countries <- c("ita", "swe", "uk", "us")
genders <- c("man", "woman")
methods <- c("exact", "euler", "3ord")
lambdas <- c(NA, 1)

grid.values <- expand.grid(countries, genders, methods, lambdas) %>% 
  setNames(c("country", "gender", "method.ssm", "lambda")) %>% 
  as.matrix()
grid.values[grid.values[, 'lambda'] == ' 1', 'lambda'] <- 1


# Rolling forecast
start <- Sys.time()
if (.Platform$OS.type == "unix"){
  
  cat("Unix OS detected. Parall execution of forecast in the 48 settings. \n")
  
  registerDoParallel(cores = min(rep,
                                 maxcl))
  cl_forecast <- makeCluster(24, type = "FORK")
  
  res_forward <- parApply(cl_forecast,
                          grid.values,
                          1,
                          function(x) 
                            rolling_uq(cg = paste(x['country'], x['gender'], sep = '_'), 
                                       method.ssm = x['method.ssm'],
                                       lambda = as.numeric(x['lambda']),
                                       n_for = n_for,
                                       parallel = FALSE))
  
} else {
  cat("Windows OS detected. Sequential execution of forecast in the 48 settings. \n")
  
  res_forward <- apply(grid.values,
                       1,
                       function(x) 
                         rolling_uq(cg = paste(x['country'], x['gender'], sep = '_'), 
                                    method.ssm = x['method.ssm'],
                                    lambda = as.numeric(x['lambda']),
                                    n_for = n_for,
                                    parallel = FALSE))
}
end <- Sys.time()


names(res_forward) <- apply(grid.values, 1, paste, collapse = "_")

output_collector <- append(output_collector, 
                           list(forecast_settings = grid.values,
                                execution_time = end - start,
                                warnings = warnings()))

save(list = c('output_collector',
              'res_forward'),
     file = here('output', '6models_forecast.Rdata'))
