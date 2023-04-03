################################################################
################################################################
## R-code for estimating CP-splines as presented in the paper
## "Smooth Constrained Mortality Forecasting"
## ® by Carlo G. Camarda
## 2019.08.28
## published in Demographic Research 

## example prepared on US males mortality data 
## ages 0-105 observed from 1960 to 2016 and forecast up to 2050

## small changes allow to model and forecast all populations presented in the manuscript

## data source: Human Mortality Database

## R-version used: R 3.6.1

## required package: 
## - MortalitySmooth
## - HMDHFDplus 

## required external files: 
## - SmoothConstrainedMortalityForecasting_Functions.R
## - SmoothConstrainedMortalityForecasting_LifeTableFunctions.R

## .R files with functions should be stored in the same folder of the 
## current R-file, otherwise path needs to be changed

## The following program follows the procedure as described
## in the associated manuscript, and itemized in Section 3.2.2

## Additionally, in order to visualize the outcomes,
## we provided some snippet for producing 
## (really simple) plots

################################################################
################################################################

## (eventually) clean the workspace from all objects
rm(list = ls())

## loading useful packages
library(MortalitySmooth) 
## for the following routines:
##  MortSmooth_tpower()
##  MortSmooth_BWB()
##  MortSmooth_BcoefB()
##  cleversearch() borrowed from library(svcm)
library(HMDHFDplus)
## for extracting data directly from the website of the Human Mortality Database by
##  readHMDweb()

## loading functions associated with CP-splines
source("SmoothConstrainedMortalityForecasting_Functions.R")
## loading functions associated with lifetable construction
source("SmoothConstrainedMortalityForecasting_LifeTableFunctions.R")

## ages
ages <- 0:105
m <- length(ages)
## observed years
years1 <- 1960:2016
n1 <- length(years1)
## all years: observed+future
years <- years1[1]:2050
n <- length(years)

## select population
pop <- "USA"
## !! Please change the object "pop"
##    if you aim to model a different population 
##    from the HMD. 
##    We tested for all populations presented in the paper:
##    - United States: USA
##    - Denmark: DNK
##    - Japan: JPN
##    - France: FRACNP

## select sex
sex <- "Male" ## altermative "Female"
## We tested both sexes for the previous populations


## loading exposures from HMD 
Ehmd <- readHMDweb(CNTRY = pop, 
                   item = "Exposures_1x1", 
                   fixup = TRUE)
## loading deaths from HMD 
Yhmd <- readHMDweb(CNTRY = pop, 
                   item = "Deaths_1x1",
                   fixup = TRUE)
## in both two previous steps
## if arguments "username" and "password" are not
## provided you will be prompted about  

## selecting ages and years
E0 <- subset(Ehmd, Year%in%years1 & Age%in%ages, sex)
## place in a mXn1 matrix
E1 <- matrix(E0[,1], m, n1)

## selecting ages and years
Y0 <- subset(Yhmd, Year%in%years1 & Age%in%ages, sex)
## place in a mXn1 matrix
Y1 <- matrix(Y0[,1], m, n1)

## observed log-rates
ETA1 <- log(Y1/E1)

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
  FITinf$bic
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
                      verbose=TRUE)
## estimated and forecast linear predictor, log-mortality
ETA.hatC <- FITcon$ETA

## extract deviance residuals
res <- FITinf$res
## replace res == NA (where Y1==0) with random residuals
## anyway they will be taken randomly in the next step
whi0 <- which(is.na(res))
res1 <- c(res)[-whi0]
res[whi0] <- sample(res1, size=length(whi0))
## small correction of actual deaths when equal to zero
## to avoid computational issues when log is taken and we numerically find the root of the function
Y11 <- Y1
Y11[whi0] <- 1e-6


## bootstrapping procedure 

## in this illustrative program, for shortening the computational time,
## we set only 10 instances in the boostrapping 

## NB: each step takes about 5 seconds on a 
##     portable personal computer, 
##     Intel i5-6300U processor, 
##     2.4 GHz $\times$ 4 and 16 Gbytes random-access memory


## number of instances
n.sim <- 10
## empty array for bootstrap deaths
Y1s <- array(0, dim=c(m,n1,n.sim))
## function to invert the bootstrap-deaths
InvDev <- function(X, Z, C){X - Z * log(X) - C}

## resample deaths
for(k in 1:n.sim){
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
for(k in 1:n.sim){
  Y1.k <- Y1s[,,k]
  FITinf.k <- PSinfant(Y=Y1.k, E=E1, lambdas=OPTinf$par, WEI=WEI1)
  deltas.k <- deltasFUN(FITinf.k)
  DELTAs[[k]] <- deltas.k
  cat(k, "\n")
}


## fit and forecast each sampled/boostrap matrix of death counts
## with CP-splines
ETA.hatCs <- array(0, dim=c(m, n, n.sim))
for(k in 1:n.sim){
  Y.k <- Y
  Y.k[,1:n1] <- Y1s[,,k]
  FITcon.k <- CPSfunction(Y=Y.k, E=E, lambdas=OPTinf$par,
                          WEI=WEI, deltas=DELTAs[[k]], S=S,
                          st.alphas=FITcon$ALPHAS)
  ETA.hatCs[,,k] <- FITcon.k$ETA
  cat(k, FITcon.k$it, FITcon.k$dalphas, "\n")
}


## computing life expectancy and lifespan measure

## observed mortality
e0.act <- numeric(n1)
ed.act <- numeric(n1)
for(i in 1:n1){
  LT <- lifetable(ages, Nx=E1[,i], Dx=Y1[,i], sex=strsplit(sex,"")[[1]][1])
  e0.act[i] <- LT$ex[1]
  ed.act[i] <- eDagger(LT)[1]
}
## estimated and forecast mortality for each boostrap instance
e0.hatC <- matrix(0, n, n.sim)
ed.hatC <- matrix(0, n, n.sim)
for(i in 1:n){
  for(k in 1:n.sim){
    mx.ik <- exp(ETA.hatCs[,i,k])
    LT <- lifetable.mx(ages, mx=mx.ik, sex=strsplit(sex,"")[[1]][1])
    e0.hatC[i,k] <- LT$ex[1]
    ed.hatC[i,k] <- eDagger(LT)[1]
  }
}




## some simple plot
## some simple plot
## some simple plot
## some simple plot
## some simple plot



## plotting relative derivatives

## smooth age-patterns and rate-of-aging
par(mfrow=c(1,2))
matplot(ages, ETA1.hatI, t="l", col=8, lty=1,
        xlab="age", ylab="log-mortality")
eta1.a.up <- apply(ETA1.hatI, 1, quantile, probs=0.975)
eta1.a.low <- apply(ETA1.hatI, 1, quantile, probs=0.025)
lines(ages, eta1.a.up, col=2, lwd=2)
lines(ages, eta1.a.low, col=2, lwd=2)
matplot(ages, FITinf$ETA1a, t="l", col=8, lty=1,
        xlab="age", ylab="rate-of-aging")
lines(ages, deltas$delta.a.up, lwd=2, col=2)
lines(ages, deltas$delta.a.low, col=2, lwd=2)


## smooth time-trends and rate-of-change
par(mfrow=c(1,2))
matplot(years1, t(ETA1.hatI), t="l", col=8, lty=1,
        xlab="year", ylab="log-mortality")

matplot(ages, FITinf$ETA1t, t="l", col=8, lty=1,
        xlab="age", ylab="rate-of-change")
lines(ages, deltas$delta.t.up, lwd=2, col=2)
lines(ages, deltas$delta.t.low, col=2, lwd=2)

## take two ages (50, 80) as example
par(mfrow=c(2,2))
plot(years1, ETA1[51,], main="age=50",
     xlab="year", ylab="log-mortality")
lines(years1, ETA1.hatI[51,], col=2, lwd=2)
plot(years1, FITinf$ETA1t[51,], t="l", col=2, lwd=2,
     xlab="year", ylab="rate-of-change")
abline(h=c(deltas$delta.t.up[51], deltas$delta.t.low[51]), col=4, lty=2)
abline(h=0, col=8, lty=3)
plot(years1, ETA1[81,], main="age=80",
     xlab="year", ylab="log-mortality")
lines(years1, ETA1.hatI[81,], col=2, lwd=2)
plot(years1, FITinf$ETA1t[81,], t="l", col=2, lwd=2,
     xlab="year", ylab="rate-of-change")
abline(h=c(deltas$delta.t.up[81], deltas$delta.t.low[81]), col=4, lty=2)
abline(h=0, col=8, lty=3)

## actual, estimated and forecast log-mortality
 
## selected ages over years
whi <- floor(seq(1,m,length=4))
par(mfrow=c(2,2))
for(i in 1:4){
  ## actual
  plot(years1, ETA1[whi[i],], xlim=range(years), 
       ylim=range(ETA1[whi[i],], ETA.hatC[whi[i],], na.rm=TRUE, finite=TRUE),
       xlab="year", ylab="log-mortality",
       main=paste("age:", ages[whi[i]]))
  ## all n.sim boostrap instances
  matlines(years, ETA.hatCs[whi[i],,], col=8, lty=1)
  ## estimated values
  lines(years, ETA.hatC[whi[i],], lwd=2, col=2)
}


## selected years over ages
whi <- floor(seq(1,n,length=4))
rany <- c(-11, 1)
par(mfrow=c(2,2))
for(i in 1:4){
  ## actual
  plot(ages, ETA1[,1], ylim=rany, t="n",
       xlab="age", ylab="log-mortality",
       main=paste("year:", years[whi[i]]))
  if(whi[i]<=n1){
    points(ages, ETA1[,whi[i]])
  }
  ## all n.sim boostrap instances
  matlines(ages, ETA.hatCs[,whi[i],], col=8, lty=1)
  ## estimated values
  lines(ages, ETA.hatC[, whi[i]], lwd=2, col=2)
}


## plotting summary measures
par(mfrow=c(1,2))
plot(years1, e0.act, xlim=range(years), ylim=range(e0.act, e0.hatC),
     xlab="year", ylab="life expactancy")
matlines(years, e0.hatC, t="l", lty=1, col=8)

plot(years1, ed.act, xlim=range(years), ylim=range(ed.act, ed.hatC, na.rm=TRUE),
     xlab="year", ylab="e-dagger")
matlines(years, ed.hatC, t="l", lty=1, col=8)
