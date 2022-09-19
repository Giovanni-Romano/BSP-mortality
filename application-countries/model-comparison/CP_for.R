require(here)
require(tidyverse)
require(MortalitySmooth) 
require(doParallel)

rm(list=ls())

set.seed(3415)

load(here('output','mortality.Rdata'))
output_collector <- list()
output_collector <- append(output_collector, list(info=Sys.time()))

# Loading function needed to fit the model 
# (avaialble from the attached material to the original paper)
source(here('application-countries',
            'model-comparison',
            'CP',
            'SmoothConstrainedMortalityForecasting_Functions.R'))
source(here('application-countries',
            'model-comparison',
            'CP',
            'SmoothConstrainedMortalityForecasting_LifeTableFunctions.R'))

# Set up minimum training window, how many step-ahead to predict
train <- 58
h_step <- 10
output_collector <- append(output_collector, list(train = train, h_step = h_step))

years0 <- years

# Function to fit  and forecast
fit_sim <- function(t, train, Y0, E0, years0, nsim = 10){
  m <- length(ages)
  years1 <- years0[1:t]
  n1 <- length(years1)
  ## all years: observed+future
  years <- years0[1:(t+h_step)]
  n <- length(years)
  
  Y1 <- Y0[,1:n1]
  E1 <- E0[,1:n1]
  
  ## observed log-rates
  ETA1 <- log(Y1/E1) #rows age, cols years
  
  ## arbitrary values for forecasting
  Y <- matrix(10, m, n)
  Y[1:m,1:n1] <- Y1
  E <- matrix(10, m, n)
  E[1:m,1:n1] <- E1
  ## 0/1 weights for the arbitray values, V in the manuscript
  WEI1 <- matrix(1, m, n1)
  ## place zero when expsoures are equal to zero
  WEI1[E1==0] <- 0
  WEI <- cbind(WEI1, matrix(0, m, n-n1))
  ## function for simply extract BIC, given deaths, exposures and weights
  BICinf <- function(par){
    FITinf <- PSinfant(Y=Y1, E=E1, lambdas=par, WEI=WEI1)
    return(FITinf$bic)
  }
  ## optimizing lambdas using greedy grid search
  OPTinf <- cleversearch(BICinf, lower=c(-4, 1), upper=c(0, 5),
                         ngrid=5, logscale=TRUE, verbose=FALSE)
  ## estimating mortality with optimal lambdas
  FITinf <- PSinfant(Y=Y1, E=E1, lambdas=OPTinf$par, WEI=WEI1)
  ## extract estimated linear predictor, log-mortality
  ETA1.hatI <- FITinf$ETA
  
  ## extract deltas from PSinfant() fitted object 
  deltas <- deltasFUN(FITinf)
  
  ## where to apply the constraints
  S <- matrix(1, m, n)
  S[,1:n1] <- 0
  
  ## modelling with CP-splines
  FITcon <- CPSfunction(Y=Y, E=E, lambdas=OPTinf$par,
                        WEI=WEI, deltas=deltas, S=S,
                        verbose=FALSE)
  ## estimated and forecast linear predictor, log-mortality
  ETA.hatC <- FITcon$ETA
  
  rates_pred <- tail(x = ETA.hatC, n = c(m,h_step))
  colnames(rates_pred) <- years0[(t+1):(t+h_step)]
  
  return(rates_pred)
}

# Function to set up the rolling window: minimum window given by train, then 
# increasing the training set by 1 each time
rolling <- function(cg){
  country <- sub("_.*", "", cg)
  gender <- sub(".*_", "", cg)
  print(paste('Doing', country, gender))
  Y0 <- eval(parse(text = paste('Y', country, gender, sep = '_')))
  E0 <- eval(parse(text = paste('N', country, gender, sep = '_')))
  
  ## observed years
  years0 <- years0[1:nrow(Y0)]

  CPsim <- lapply(train:length(years0),
                  . %>%
                    fit_sim(., 
                            train = train, 
                            Y0 = t(Y0),
                            E0 = t(E0),
                            years0 = years0,
                            nsim = 100))
  
  return(CPsim)
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
     file = here('output','CP_for.Rdata'))


