require(here)
require(magrittr)
require(purrr)
require(dplyr)
require(MortalitySmooth) 
require(parallel)

rm(list=ls())

set.seed(3415)
RNGkind("L'Ecuyer-CMRG")

options(mc.cores = 10)

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
fit_sim <- function(t, train, Y0, E0, years0, nsim){
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
  
  ## Bootstrap
  res <- FITinf$res
  whi0 <- which(is.na(res))
  res1 <- c(res)[-whi0]
  res[whi0] <- sample(res1, size=length(whi0))
  Y11 <- Y1
  Y11[whi0] <- 1e-6
  Y1s <- array(0, dim=c(m,n1,nsim))
  ## function to invert the bootstrap-deaths
  InvDev <- function(X, Z, C){X - Z * log(X) - C}
  ## resample deaths
  for(k in 1:nsim){
    ## resample (with replacement) the residuals
    resS <- matrix(sample(x=res, size=m*n1, replace = TRUE), m, n1)  
    ## computing C =  res^2/2 + Dth - Dth*ln(Dth)
    C <- 1/2*(resS^2) + Y11 - Y11*log(Y11)
    ## recomputing the deaths
    LOW <- Y11
    LOW[which(resS>0)] <- 0
    UP <- Y11
    UP[which(resS<0)] <- 1e8
    UP <- ifelse(LOW>=UP, UP+10^-8, UP)
    Y1s.i <- matrix(NA,m,n1)
    for(i in 1:m){
      for(j in 1:n1){
        Y1s.i[i,j]  <- uniroot(InvDev, lower=LOW[i,j], upper=UP[i,j],
                               tol = 0.00001, Z=Y11[i,j], C=C[i,j])$root
      }
    }
    Y1s[,,k] <- Y1s.i
    cat(k, "\n")
  }
  ## fit each sampled/bootstrap matrix of death counts
  ## and from the derivatives, compute associated confidence intervals
  DELTAs <- list()
  for(k in 1:nsim){
    Y1.k <- Y1s[,,k]
    FITinf.k <- PSinfant(Y=Y1.k, E=E1, lambdas=OPTinf$par, WEI=WEI1)
    deltas.k <- deltasFUN(FITinf.k)
    DELTAs[[k]] <- deltas.k
    cat(k, "\n")
  }
  ## fit and forecast each sampled/boostrap matrix of death counts
  ## with CP-splines
  ETA.hatCs <- array(0, dim=c(m, n, nsim))
  for(k in 1:nsim){
    Y.k <- Y
    Y.k[,1:n1] <- Y1s[,,k]
    FITcon.k <- CPSfunction(Y=Y.k, E=E, lambdas=OPTinf$par,
                            WEI=WEI, deltas=DELTAs[[k]], S=S,
                            st.alphas=FITcon$ALPHAS)
    ETA.hatCs[,,k] <- FITcon.k$ETA
    cat(k, FITcon.k$it, FITcon.k$dalphas, "\n")
  }
  ##
  
  ETA.hatCs_for <- tail(x = ETA.hatCs, n = c(m,h_step,nsim))
  ETA.hatC_U95 <- apply(ETA.hatCs_for, MARGIN = c(1,2), function(x){quantile(x, 0.975)})
  colnames(ETA.hatC_U95) <- years0[(t+1):(t+h_step)]
  ETA.hatC_L95 <- apply(ETA.hatCs_for, MARGIN = c(1,2), function(x){quantile(x, 0.025)})
  colnames(ETA.hatC_L95) <- years0[(t+1):(t+h_step)]
  
  rates_pred <- tail(x = ETA.hatC, n = c(m,h_step))
  colnames(rates_pred) <- years0[(t+1):(t+h_step)]
  
  # return(rates_pred)
  return(list(pred = rates_pred,
              U95 = ETA.hatC_U95,
              L95 = ETA.hatC_L95))
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
  
  CPsim <- mclapply(train:length(years0),
                    FUN = . %>%
                      fit_sim(., 
                              train = train, 
                              Y0 = t(Y0),
                              E0 = t(E0),
                              years0 = years0,
                              nsim = 500))
  
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
     file = here('output', 'CP_boot_for.Rdata'))


