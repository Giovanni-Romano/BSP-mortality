require(here)
require(tidyverse)
require(KFAS)
require(doParallel)

rm(list=ls())

set.seed(4238)

source(here('source','BSP.R'))
source(here('source','setup.R'))
load(here('output','mortality.Rdata'))

# Data preparation
data <- Y_ita_woman/N_ita_woman

collector <- list()
collector <- append(collector, list(kernel = matern_kernel,
                                    delta = delta,
                                    age_knots = age_knots,
                                    rep = rep))

params <- expand.grid(phi = c(0.1,0.5,1,1.5), 
                      kappa = c(1,1.5,2,2.5))
res <- matrix(NA, nrow = nrow(params), ncol = 4)
for(it in 1:nrow(params)){
  phi <- params$phi[it]
  kappa <- params$kappa[it]
  kernel <- function(x){
    matern_kernel(x, phi = phi, kappa = kappa)
  }
  model <- bsp.model(rates = data,
                     delta = delta,
                     age_knots = age_knots,
                     kernel = kernel)
  fit <- bsp.fit(model, 
                 rep = 1, 
                 method = 'Nelder-Mead', 
                 parallel = FALSE, 
                 maxcl = 50)
  res[it,] <- fit$info$optim$par
}

save.image(here("output", "robustness.Rdata"))
