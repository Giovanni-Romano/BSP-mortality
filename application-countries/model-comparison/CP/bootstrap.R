## bootstrapping procedure 

## in this illustrative program, for shortening the computational time,
## we set only 10 instances in the boostrapping 

## NB: each step takes about 5 seconds on a 
##     portable personal computer, 
##     Intel i5-6300U processor, 
##     2.4 GHz $\times$ 4 and 16 Gbytes random-access memory


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

## number of instances
n.sim <- 100
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
