################################################################
################################################################
## R-code with a set of function useful 
## for building life-table based on different inputs
## and extract e-dagger from it

## Associated to the paper (and Supplementary Material)
## "Smooth Constrained Mortality Forecasting"
## ® by Carlo G. Camarda
## 2019.08.28
## published in Demographic Research
## Volume 41, Article 38, pages 1091-1130. 
## Published 24 October 2019. 
## http://www.demographic-research.org/Volumes/Vol41/38/
## DOI: 10.4054/DemRes.2019.41.38

## function for constructing a classic (& rather general) lifetable
## given a set of death and exposures
lifetable <- function(x, Nx, Dx, sex="M", ax=NULL){
  m <- length(x)
  mx  <- Dx/Nx
  n <- c(diff(x), NA)
  if(is.null(ax)){
    ax <- rep(0,m)
    if(x[1]!=0 | x[2]!=1){
      ax <- n/2
      ax[m] <- 1 / mx[m]
    }else{    
      if(sex=="F"){
        if(mx[1]>=0.107){
          ax[1] <- 0.350
        }else{
          ax[1] <- 0.053 + 2.800*mx[1]
        }
      }
      if(sex=="M"){
        if(mx[1]>=0.107){
          ax[1] <- 0.330
        }else{
          ax[1] <- 0.045 + 2.684*mx[1]
        }
      }
      ax[-1] <- n[-1]/2
      ax[m] <- 1 / mx[m]
    }
  }
  qx  <- n*mx / (1 + (n - ax) * mx)
  qx[m] <- 1
  px  <- 1-qx
  lx  <- cumprod(c(1,px))*100000
  dx  <- -diff(lx)
  Lx  <- n*lx[-1] + ax*dx
  lx <- lx[-(m+1)]
  Lx[m] <- lx[m]/mx[m]
  Lx[is.na(Lx)] <- 0 ## in case of NA values
  Lx[is.infinite(Lx)] <- 0 ## in case of Inf values
  Tx  <- rev(cumsum(rev(Lx)))
  ex  <- Tx/lx
  return.df <- data.frame(x, n, Nx, Dx, mx, ax, qx, px, lx, dx, Lx, Tx, ex)
  return(return.df)
}

## function for constructing a lifetable starting from rates
lifetable.mx <- function(x, mx, sex="M", ax=NULL){
  m <- length(x)
  n <- c(diff(x), NA)
  if(is.null(ax)){
    ax <- rep(0,m)
    if(x[1]!=0 | x[2]!=1){
      ax <- n/2
      ax[m] <- 1 / mx[m]
    }else{    
      if(sex=="F"){
        if(mx[1]>=0.107){
          ax[1] <- 0.350
        }else{
          ax[1] <- 0.053 + 2.800*mx[1]
        }
      }
      if(sex=="M"){
        if(mx[1]>=0.107){
          ax[1] <- 0.330
        }else{
          ax[1] <- 0.045 + 2.684*mx[1]
        }
      }
      ax[-1] <- n[-1]/2
      ax[m] <- 1 / mx[m]
    }
  }
  qx  <- n*mx / (1 + (n - ax) * mx)
  qx[m] <- 1
  px  <- 1-qx
  lx  <- cumprod(c(1,px))*100000
  dx  <- -diff(lx)
  Lx  <- n*lx[-1] + ax*dx
  lx <- lx[-(m+1)]
  Lx[m] <- lx[m]/mx[m]
  Lx[is.na(Lx)] <- 0 ## in case of NA values
  Lx[is.infinite(Lx)] <- 0 ## in case of Inf values
  Tx  <- rev(cumsum(rev(Lx)))
  ex  <- Tx/lx
  return.df <- data.frame(x, n, mx, ax, qx, px, lx, dx, Lx, Tx, ex)
  return(return.df)
}

## function for computing e-dagger from a lifetable object
## obtained from the previous routines
eDagger <- function(LT){
  m <- nrow(LT)
  m1 <- m-1
  ## lost length
  LT$lost <- c(LT$ex[-1]+1-LT$ax[-m], NA)
  ## positive lx
  LT$lxpos <- ifelse(abs(LT$lx)<0.001, 0.05, LT$lx)
  ## compute e-dagger for each x
  ed <- numeric(m)
  for(i in 1:m){
    p1 <- sum(LT$dx[i:m1] * LT$lost[i:m1], na.rm=TRUE)
    p2 <- 0.5*LT$ex[m]*LT$lx[m]
    p3 <- LT$lxpos[i]
    ed[i] <- (p1+p2)/p3
  }
  return(ed)
}
