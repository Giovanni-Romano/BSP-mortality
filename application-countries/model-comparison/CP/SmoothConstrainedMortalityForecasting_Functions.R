################################################################
################################################################
## R-code with a set of functions useful 
## for modelling and forecasting 
## based on CP-splines

## Associated to the paper (and Supplementary Material)
## "Smooth Constrained Mortality Forecasting"
## ® by Carlo G. Camarda
## 2019.08.28
## published in Demographic Research
## Volume 41, Article 38, pages 1091-1130. 
## Published 24 October 2019. 
## http://www.demographic-research.org/Volumes/Vol41/38/
## DOI: 10.4054/DemRes.2019.41.38

## function to build up B-splines and associated bases for derivatives
BsplineGrad <- function(x, xl, xr, ndx=NULL, deg, knots=NULL){
  if(is.null(knots)){
    dx <- (xr - xl)/ndx
    knots <- seq(xl - deg * dx, xr + deg * dx, by=dx)
    knots <- round(knots, 8)
  }else{
    knots <- knots
    dx <- diff(knots)[1]
  }
  P <- outer(x, knots, MortSmooth_tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff=deg+1)/(gamma(deg+1)*dx^deg)
  B <- (-1)^(deg + 1) * P %*% t(D)
  ##
  knots1 <- knots[-c(1,length(knots))]
  P <- outer(x, knots1, MortSmooth_tpower, deg-1)
  n <- dim(P)[2]
  D <- diff(diag(n),diff=deg)/(gamma(deg)*dx^(deg-1))
  BB <- ((-1)^(deg) * P %*% t(D))/dx
  D <- diff(diag(ncol(BB) + 1))
  C <- BB %*% D
  ##
  out <- list(dx=dx, knots=knots, B=B, C=C)
}


## function for estimating two-dimensional P-spline 
## addressing infant mortality (when infant=TRUE)
## for a given set of smoothing parameters lambdas
PSinfant <- function(Y, E, lambdas, WEI, infant=TRUE, verbose=FALSE){
  ## dimensions
  m <- nrow(Y)
  n <- ncol(Y)
  a <- 1:m
  t <- 1:n
  ## w/o age 0
  a0 <- a[-1]
  m0 <- m-1
  ## original offset
  OFF <- log(E)
  ## B-splines basis
  ## with infant-specilized coeff
  if(infant){
    ## over ages w/o age 0
    a0min <- min(a0)
    a0max <- max(a0)
    nda0 <- floor(m0/5)
    dega <- 3
    BCa0 <- BsplineGrad(a0, a0min, a0max, nda0, dega)
    Ba0 <- BCa0$B
    nba0 <- ncol(Ba0)
    ## adding infant-specific basis
    Ba <- cbind(0, Ba0)
    Ba <- rbind(c(1, rep(0,nba0)), Ba)
    nba <- ncol(Ba)
    ## basis for the derivatives
    ## over ages w/o age 0
    Ca0 <- BCa0$C
    ## adding infant-specific basis
    Ca <- cbind(0, Ca0)
    Ca <- rbind(c(-1,Ba[2,1:dega+1],rep(0,nba0-dega)),
                Ca)
  }else{
    amin <- min(a)
    amax <- max(a)
    nda <- floor(m/5)
    dega <- 3
    BCa <- BsplineGrad(a, amin, amax, nda, dega)
    Ba <- BCa$B
    nba <- ncol(Ba)
    ## basis for the derivatives
    Ca <- BCa$C
  }
  ## over years
  tmin <- min(t)
  tmax <- max(t)
  ndt <- floor(n/5)
  degt <- 3
  BCt <- BsplineGrad(t, tmin, tmax, ndt, degt)
  Bt <- BCt$B
  nbt <- ncol(Bt)
  ## basis for the derivatives
  Ct <- BCt$C
  ## weights for exposures=0
  WEI[E==0] <- 0
  
  ## tensor product of B-splines for the GLAM
  Ba1 <- kronecker(matrix(1, ncol=nba, nrow=1), Ba)
  Ba2 <- kronecker(Ba, matrix(1, ncol=nba, nrow=1))
  RTBa <- Ba1 * Ba2
  Bt1 <- kronecker(matrix(1, ncol=nbt, nrow=1), Bt)
  Bt2 <- kronecker(Bt, matrix(1, ncol=nbt, nrow=1))
  RTBt <- Bt1 * Bt2

  ## penalty terms
  ## over ages
  Da <- diff(diag(nba), diff=2)
  ## no penalization over age for age 0,
  ## with infant=TRUE
  if(infant){
    Da[1,1] <- 0
  }
  tDDa <- t(Da) %*% Da
  ## over years
  Dt <- diff(diag(nbt), diff=2)
  tDDt <- t(Dt) %*% Dt
  ## kronecker product of difference matrices
  Pa <- kronecker(diag(nbt), tDDa)
  Pt <- kronecker(tDDt, diag(nba))
  ## smoothing parameters
  lambda.a <- lambdas[1]
  lambda.t <- lambdas[2]
  P <- lambda.a * Pa + lambda.t * Pt
  ## data in vector for the starting values
  y <- c(Y)
  e <- c(E)
  wei <- c(WEI)
  off <- log(e)
  off0 <- off
  off0[wei==0] <- 100
  ## "other" offset in matrices
  OFF0 <- matrix(off0, m, n)
  ## starting coeff
  aa <- rep(a, n)
  tt <- rep(t, each = m)
  fit0 <- glm(round(y) ~ aa + tt + offset(off0), 
              family = poisson, weights = wei)
  etaGLM <- matrix(log(fit0$fitted) - c(off0), m, n)
  eta0 <- log((Y + 1)) - OFF0
  eta0[WEI == 0] <- etaGLM[WEI == 0]
  BBa <- solve(t(Ba) %*% Ba + 1e-06*diag(nba), t(Ba))
  BBt <- solve(t(Bt) %*% Bt + 1e-06*diag(nbt), t(Bt))
  alphas <- MortSmooth_BcoefB(BBa, BBt, eta0)
  ## Poisson iteration within a GLAM framework
  for(it in 1:20){
    eta <- MortSmooth_BcoefB(Ba, Bt, alphas)
    mu <- exp(OFF0 + eta)
    W <- mu
    z <- eta + (1/mu) * (Y - mu)
    z[which(WEI == 0)] <- 0
    WW <- WEI * W
    BWB <- MortSmooth_BWB(RTBa, RTBt, nba, nbt, WW)
    BWBpP <- BWB + P
    BWz <- MortSmooth_BcoefB(t(Ba), t(Bt), (WW * z))
    alphas0 <- solve(BWBpP, c(BWz))
    alphas.old <- alphas
    alphas <- matrix(alphas0, nrow = nba)
    dalphas <- max(abs(alphas-alphas.old)/abs(alphas))
    if(verbose) cat(it, dalphas, "\n")
    if(dalphas<=10^-6 & it>=4) break
  }
  ## fitted values
  ALPHAS.hat <- alphas
  ETA.hat <- eta
  Y.hat <- exp(OFF0 + ETA.hat)
  ## fitted derivatives over ages
  ETA.hat1a <- MortSmooth_BcoefB(Ca, Bt, ALPHAS.hat)
  ## fitted derivatives over years
  ETA.hat1t <- MortSmooth_BcoefB(Ba, Ct, ALPHAS.hat)
  ## diagnostics
  H <- solve(BWBpP, BWB)
  h <- diag(H)
  ed <- sum(h)
  y0 <- y
  y0[y==0] <- 10^-8
  dev <- 2 * sum(wei * (y * log(y0/mu)))
  bic <- dev + log(sum(wei))*ed
  ## deviance residuals
  res0 <- sign(Y - Y.hat)
  res1 <- log(Y/Y.hat)
  res2 <- Y - Y.hat
  res <- res0 * sqrt(2 * (Y * res1 - res2))
  res[which(is.nan(res))] <- NA
  res <- matrix(res,m,n)

  ## return objects
  out <- list(
    ## original data
    Y=Y, E=E, lambdas=lambdas,
    ## diagnostic
    dev=dev, ed=ed, bic=bic, res=res,
    ## fitted values
    ALPHAS=ALPHAS.hat, ETA=ETA.hat, MU=Y.hat,
    ETA1a=ETA.hat1a, ETA1t=ETA.hat1t, 
    ## basis and associated objects
    BCa=ifelse(infant, BCa0, BCa), BCt=BCt,
    ## convergence objects
    dalphas=dalphas, it=it
  )
  return(out)
}


## function to extract deltas from a PSinfant object
## for a given confidence level, default 95% and 50% 
## over age and year as reccommened/used in the manuscript
deltasFUN <- function(object, levels=c(95,50)){
  ETA1a <- object$ETA1a
  ETA1t <- object$ETA1t
  ## compute levels and deltas over ages
  p.a.up <- (100 - (100-levels[1])/2)/100
  p.a.low <- ((100-levels[1])/2)/100
  delta.a.up <- apply(ETA1a, 1, quantile,
                      probs=p.a.up)
  delta.a.low <- apply(ETA1a, 1, quantile,
                       probs=p.a.low)
  ## compute levels and deltas over years
  p.t.up <- (100 - (100-levels[2])/2)/100
  p.t.low <- ((100-levels[2])/2)/100
  delta.t.up <- apply(ETA1t, 1, quantile,
                      probs=p.t.up)
  delta.t.low <- apply(ETA1t, 1, quantile,
                       probs=p.t.low)
  ## return objects
  out <- list(delta.a.up=delta.a.up,
              delta.a.low=delta.a.low,
              delta.t.up=delta.t.up,
              delta.t.low=delta.t.low)
  return(out)
}


## function for estimating CP-splines
## for a given set of smoothing parameters lambdas
## commonly taken from the estimated object, PSinfant()
CPSfunction <- function(Y, E, lambdas, WEI,
                        kappas=c(10^4, 10^4),
                        deltas, S,
                        infant=TRUE, verbose=FALSE,
                        st.alphas=NULL){
  ## dimensions
  m <- nrow(Y)
  n <- ncol(Y)
  a <- 1:m
  t <- 1:n
  ## w/o age 0
  a0 <- a[-1]
  m0 <- m-1
  ## original offset
  OFF <- log(E)
  ## B-splines basis
  ## with infant-specilized coeff
  if(infant){
    ## over ages w/o age 0
    a0min <- min(a0)
    a0max <- max(a0)
    nda0 <- floor(m0/5)
    dega <- 3
    BCa0 <- BsplineGrad(a0, a0min, a0max, nda0, dega)
    Ba0 <- BCa0$B
    nba0 <- ncol(Ba0)
    ## adding infant-specific basis
    Ba <- cbind(0, Ba0)
    Ba <- rbind(c(1, rep(0,nba0)), Ba)
    nba <- ncol(Ba)
    ## basis for the derivatives
    ## over ages w/o age 0
    Ca0 <- BCa0$C
    ## adding infant-specific basis
    Ca <- cbind(0, Ca0)
    Ca <- rbind(c(-1,Ba[2,1:dega+1],rep(0,nba0-dega)),
                Ca)
  }else{
    amin <- min(a)
    amax <- max(a)
    nda <- floor(m/5)
    dega <- 3
    BCa <- BsplineGrad(a, amin, amax, nda, dega)
    Ba <- BCa$B
    nba <- ncol(Ba)
    ## basis for the derivatives
    Ca <- BCa$C
  }
  ## over years
  tmin <- min(t)
  tmax <- max(t)
  ndt <- floor(n/5)
  degt <- 3
  BCt <- BsplineGrad(t, tmin, tmax, ndt, degt)
  Bt <- BCt$B
  nbt <- ncol(Bt)
  ## basis for the derivatives
  Ct <- BCt$C
  ## weights for exposures=0
  WEI[E==0] <- 0
  
  ## tensor product of B-splines for the GLAM
  Ba1 <- kronecker(matrix(1, ncol=nba, nrow=1), Ba)
  Ba2 <- kronecker(Ba, matrix(1, ncol=nba, nrow=1))
  RTBa <- Ba1 * Ba2
  Bt1 <- kronecker(matrix(1, ncol=nbt, nrow=1), Bt)
  Bt2 <- kronecker(Bt, matrix(1, ncol=nbt, nrow=1))
  RTBt <- Bt1 * Bt2
  ## tensor product for the derivatives
  Ca1 <- kronecker(matrix(1, ncol=nba, nrow=1), Ca)
  Ca2 <- kronecker(Ca, matrix(1, ncol=nba, nrow=1))
  RTCa <- Ca1 * Ca2
  Ct1 <- kronecker(matrix(1, ncol=nbt, nrow=1), Ct)
  Ct2 <- kronecker(Ct, matrix(1, ncol=nbt, nrow=1))
  RTCt <- Ct1 * Ct2

  
  ## penalty terms for smoothing
  ## over ages
  Da <- diff(diag(nba), diff=2)
  ## no penalization over age for age 0,
  ## with infant=TRUE
  if(infant){
    Da[1,1] <- 0
  }
  tDDa <- t(Da) %*% Da
  ## over years
  Dt <- diff(diag(nbt), diff=2)
  tDDt <- t(Dt) %*% Dt
  ## kronecker product of difference matrices
  Pa <- kronecker(diag(nbt), tDDa)
  Pt <- kronecker(tDDt, diag(nba))
  ## smoothing parameters
  lambda.a <- lambdas[1]
  lambda.t <- lambdas[2]
  P <- lambda.a * Pa + lambda.t * Pt

  ## penalty terms for constraints
  ## extract deltas
  delta.a.up <- deltas$delta.a.up
  delta.a.low <- deltas$delta.a.low
  delta.t.up <- deltas$delta.t.up
  delta.t.low <- deltas$delta.t.low
  ## construct G matrices
  ones12 <- matrix(1, n, 1)
  g.a.up <- kronecker(ones12, delta.a.up)
  G.a.up <- matrix(g.a.up, m, n)
  g.a.low <- kronecker(ones12, delta.a.low)
  G.a.low <- matrix(g.a.low, m, n)
  g.t.up <- kronecker(ones12, delta.t.up)
  G.t.up <- matrix(g.t.up, m, n)
  g.t.low <- kronecker(ones12, delta.t.low)
  G.t.low <- matrix(g.t.low, m, n)

  ## data in vector for the starting values
  y <- c(Y)
  e <- c(E)
  wei <- c(WEI)
  off <- log(e)
  off0 <- off
  off0[wei==0] <- 100
  ## "other" offset in matrices
  OFF0 <- matrix(off0, m, n)
  if(is.null(st.alphas)){
    ## starting coeff
    aa <- rep(a, n)
    tt <- rep(t, each = m)
    fit0 <- glm(round(y) ~ aa + tt + offset(off0), 
                family = poisson, weights = wei)
    etaGLM <- matrix(log(fit0$fitted) - c(off0), m, n)
    eta0 <- log((Y + 1)) - OFF0
    eta0[WEI == 0] <- etaGLM[WEI == 0]
    BBa <- solve(t(Ba) %*% Ba + 1e-06*diag(nba), t(Ba))
    BBt <- solve(t(Bt) %*% Bt + 1e-06*diag(nbt), t(Bt))
    alphas <- MortSmooth_BcoefB(BBa, BBt, eta0)
  }else{
    alphas=st.alphas
  }
  ## Poisson iteration with asymmetric penalty in a GLAM framework
  for(it in 1:100){
    ## over ages
    CaBt.alphas <- MortSmooth_BcoefB(Ca, Bt, alphas)
    ## up
    v.a.up <- CaBt.alphas > G.a.up
    v.a.up <- v.a.up * S
    P.a.up <- MortSmooth_BWB(RTCa, RTBt,
                            nba, nbt, v.a.up)
    P.a.up <- kappas[1] * P.a.up
    p.a.up <- MortSmooth_BcoefB(t(Ca), t(Bt),
                                v.a.up*G.a.up)
    p.a.up <- kappas[1] * p.a.up
    ## low
    v.a.low <- CaBt.alphas < G.a.low
    v.a.low <- v.a.low * S
    P.a.low <- MortSmooth_BWB(RTCa, RTBt,
                            nba, nbt, v.a.low)
    P.a.low <- kappas[1] * P.a.low
    p.a.low <- MortSmooth_BcoefB(t(Ca), t(Bt),
                                v.a.low*G.a.low)
    p.a.low <- kappas[1] * p.a.low

    P.a <- P.a.up + P.a.low
    p.a <- p.a.up + p.a.low
    ## over years
    BaCt.alphas <- MortSmooth_BcoefB(Ba, Ct, alphas)
    ## up
    v.t.up <- BaCt.alphas > G.t.up
    v.t.up <- v.t.up * S  
    P.t.up <- MortSmooth_BWB(RTBa, RTCt, nba, nbt,
                             v.t.up)
    P.t.up <- kappas[2] * P.t.up
    p.t.up <- MortSmooth_BcoefB(t(Ba), t(Ct),
                                v.t.up*G.t.up)
    p.t.up <- kappas[2] * p.t.up
    ## low
    v.t.low <- BaCt.alphas < G.t.low
    v.t.low <- v.t.low * S  
    P.t.low <- MortSmooth_BWB(RTBa, RTCt, nba, nbt,
                             v.t.low)
    P.t.low <- kappas[2] * P.t.low
    p.t.low <- MortSmooth_BcoefB(t(Ba), t(Ct),
                                v.t.low*G.t.low)
    p.t.low <- kappas[2] * p.t.low

    P.t <- P.t.up + P.t.low
    p.t <- p.t.up + p.t.low

    eta <- MortSmooth_BcoefB(Ba, Bt, alphas)
    mu <- exp(OFF0 + eta)
    W <- mu
    z <- eta + (1/mu) * (Y - mu)
    z[which(WEI == 0)] <- 0
    WW <- WEI * W
    BWB <- MortSmooth_BWB(RTBa, RTBt, nba, nbt, WW)
    BWBpP <- BWB + P + P.a + P.t
    BWz <- MortSmooth_BcoefB(t(Ba), t(Bt), (WW * z))
    BWzpP <- BWz + p.a + p.t
    alphas0 <- solve(BWBpP, c(BWzpP))
    alphas.old <- alphas
    alphas <- matrix(alphas0, nrow = nba)
    dalphas <- max(abs(alphas-alphas.old)/abs(alphas))
    if(verbose) cat(it, dalphas, "\n")
    if(dalphas<=10^-4 & it>=4) break   
  }
  ## fitted values
  ALPHAS.hat <- alphas
  ETA.hat <- eta
  Y.hat <- exp(OFF0 + ETA.hat)
  ## fitted derivatives over ages
  ETA.hat1a <- MortSmooth_BcoefB(Ca, Bt, ALPHAS.hat)
  ## fitted derivatives over years
  ETA.hat1t <- MortSmooth_BcoefB(Ba, Ct, ALPHAS.hat)
  ## diagnostics
  H <- solve(BWBpP, BWB)
  h <- diag(H)
  ed <- sum(h)
  y0 <- y
  y0[y==0] <- 10^-8
  dev <- 2 * sum(wei * (y * log(y0/mu)))
  bic <- dev + log(sum(wei))*ed

  ## return objects
  out <- list(
    ## original data
    Y=Y, E=E, lambdas=lambdas,
    kappas=kappas, deltas=deltas, S=S,
    ## diagnostic
    dev=dev, ed=ed, bic=bic,
    ## fitted values
    ALPHAS=ALPHAS.hat, ETA=ETA.hat, MU=Y.hat,
    ETA1a=ETA.hat1a, ETA1t=ETA.hat1t, 
    ## basis and associated objects
    BCa=ifelse(infant, BCa0, BCa), BCt=BCt,
    ## convergence
    dalphas=dalphas, it=it
  )
  return(out)
}



