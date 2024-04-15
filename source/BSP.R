require(magrittr)
require(KFAS)
require(splines2)
require(invgamma)
require(doParallel)


##########################################
# Utils to define bsp models via KFAS ####
##########################################
build_splines <- function(rates,
                          knots){
  
  # Maximum age
  Z <- ncol(rates)
  ages <- as.numeric(colnames(rates))
  ages_no0 <- ages[-1]
  
  # Construction of the spline basis functions
  # K = # of basis functions
  ages_degree <- 2
  S <- bSpline(ages_no0, 
               knots = age_knots, 
               degree = ages_degree,
               intercept = TRUE)
  K <- ncol(S)
  for (k in 1:K){
    S[,k] <- S[,k]/max(S[,k])
  }
  
  return(list('S' = S, 
              'K'= K,
              'Z' = Z))
}


build_Zt <- function(Z,
                     S){
  K <- ncol(S)
  
  Zt <- matrix(0, nrow = Z, ncol = 3*(K+1))
  Zt[-1, seq(4, 3*(K+1), by=3)] <- S
  Zt[1,1] <- 1
  
  return(Zt)
}


build_Tt <- function(method.ssm,
                     delta,
                     lambda,
                     K){
  # Part in common with all 3 methods
  T_tmp <- diag(1,3)
  T_tmp[1,2] <- T_tmp[2,3] <- delta*lambda
  
  # If exact or 3ord, then fill also [1, 3]
  if (method.ssm %in% c("exact", "3ord")){
    T_tmp[1,3] <- delta^2/2*lambda^2
  }
  # If euler, then [1, 3] still 0
  
  Tt <- kronecker(diag(1,K+1), T_tmp)
  
  return(Tt)
}


build_Rt <- function(method.ssm,
                     K){
  # If exact or 3ord, then R_tmp (and hence Rt) is the identity.
  # If euler, it is not.
  if (method.ssm %in% c("3ord", "exact")){
    Rt <- diag(1,3*(K+1))
  } else if (method.ssm == "euler") {
    R_tmp <- matrix(c(0, 0,
                      1, 0,
                      0, 1), 
                    nrow = 3, ncol = 2,
                    byrow = TRUE)
    Rt <- kronecker(diag(1,K+1), R_tmp) 
  }
  
  return(Rt)
}


build_Qt <- function(method.ssm,
                     sigma2_u,
                     sigma2_a, 
                     lambda,
                     delta,
                     kernel,
                     S){
  K <- ncol(S)
  ages_max <- c(0, apply(S, 2, which.max))
  
  if (method.ssm %in% c("exact", "3ord")){
    Q_tmp <- matrix(NA,3,3)
    Q_tmp[1,1] <- sigma2_u/3 * delta^3 * lambda^2 + 
      sigma2_a/20 * delta^5 * lambda^4
    Q_tmp[1,2] <- 
      Q_tmp[2,1] <- sigma2_u/2 * delta^2 * lambda + 
      sigma2_a/8 * delta^4 * lambda^3
    Q_tmp[1,3] <- 
      Q_tmp[3,1] <- sigma2_a/6 * delta^3 * lambda^2
    Q_tmp[2,2] <- sigma2_u*delta + 
      sigma2_a/3 * delta^3 * lambda^2
    Q_tmp[2,3] <- 
      Q_tmp[3,2] <- sigma2_a/2 * delta^2 * lambda
    Q_tmp[3,3] <- sigma2_a*delta
    Qt <- kronecker(diag(1,K+1), Q_tmp)
    for(I in 1:(K+1-1)){
      for(J in (I+1):(K+1)){
        rho_u <- kernel(x = abs(ages_max[J]-ages_max[I])) 
        rho_a <- 0 
        Q_tmp[1,1] <- sigma2_u/3 * rho_u * delta^3 * lambda^2 + 
          sigma2_a/20 * rho_a * delta^5 * lambda^4
        Q_tmp[1,2] <- 
          Q_tmp[2,1] <- sigma2_u/2 * rho_u*delta^2 * lambda + 
          sigma2_a/8 * rho_a*delta^4 * lambda^3
        Q_tmp[1,3] <- 
          Q_tmp[3,1] <- sigma2_a/6 * rho_a * delta^3 * lambda^2
        Q_tmp[2,2] <- sigma2_u * rho_u * delta + 
          sigma2_a/3 * rho_a*delta^3 * lambda^2
        Q_tmp[2,3] <- 
          Q_tmp[3,2] <- sigma2_a/2 * rho_a * delta^2 * lambda
        Q_tmp[3,3] <- sigma2_a*rho_a*delta
        Qt[(1+(I-1)*3):(I*3),(1+(J-1)*3):(J*3)] <- 
          Qt[(1+(J-1)*3):(J*3),(1+(I-1)*3):(I*3)] <- Q_tmp
      }
    }
  }
  
  
  if (method.ssm == "euler"){
    Q_tmp <- delta * matrix(c(sigma2_u, 0,
                              0, sigma2_a),
                            nrow = 2, ncol = 2,
                            byrow = TRUE)
    
    Qt <- kronecker(diag(1,K+1), Q_tmp)
    for(I in 1:(K+1-1)){
      for(J in (I+1):(K+1)){
        rho_u <- kernel(x = abs(ages_max[J]-ages_max[I])) 
        rho_a <- 0 
        Qt[(1+(I-1)*2):(I*2),(1+(J-1)*2):(J*2)] <- 
          Qt[(1+(J-1)*2):(J*2),(1+(I-1)*2):(I*2)] <- c(rho_u, rho_a) * Q_tmp
      }
    }
  }
  
  return(Qt)
}


build_Ht <- function(sigma2_e,
                     Z){
  # Observational variance
  Ht <- diag(sigma2_e, Z)
  
  return(Ht)
}


build_u0 <- function(rates,
                     S){
  
  # Initial values u_0 = a1 set at frequentist estimate 
  Eta <- log(rates[,-1])
  K <- ncol(S)
  U <- matrix(NA, nrow(rates), K)
  for (t in 1:nrow(rates)){
    fit <- lm(Eta[t,] ~ S + 0)
    U[t,] <- fit$coefficients
  }
  tmp_a1 <- rep(0, 3*(K+1))
  tmp_a1[1] <- log(rates[1,1])
  tmp_a1[seq(4, 3*(K+1), by=3)] <- U[1,]
  a1 <- matrix(tmp_a1, 3*(K+1), 1)
  
  #Initial variance var(u_0) = P1
  P1 <- diag(10, 3*(K+1))
  
  return(list('a1' = a1, 
              'P1' = P1))
}

######################################################
# BSP model definition via KFAS #####################
######################################################
# INPUT:
# - data: list with named element "rates" containing the mortality rates
# - delta: frequency of the observations
# - age_knots: vector of knots for splines (excluding age 0)
# - kernel: function to define correlation among spline weights
#
# RETURN: 
# - List with the model and additional info.
#
bsp.model <- function(method.ssm,
                      lambda = NA,
                      rates, 
                      delta, 
                      age_knots,
                      kernel){
  
  info <- list()
  info <- append(info, list(data.nrow = nrow(rates),
                            data.ncol = ncol(rates),
                            age_knots = age_knots,
                            kernel = kernel,
                            delta = delta))
  
  splines <- build_splines(rates = rates,
                           knots = age_knots)
  S <- splines$S
  K <- splines$K
  Z <- splines$Z
  ages_max <- c(0, apply(S, 2, which.max))
  
  
  #### State-space model
  #
  # u_t = (u_t1, ..., u_tK) stays for 
  # beta_t = (beta_t1, ..., beta_tK) in the paper.
  #
  # Model structure in KFAS:
  # log_m_t = Zt*beta_t + N(0,Ht)
  # beta_t+1 = Tt*beta_t + Rt*N(0,Qt)
  
  # NA are hyperparameters to estimate via MLE
  # If method is 3ord then we have just one variance parameter.
  #   I set sigma2_u = 0 in this case and leave sigma2_a free even if there 
  #   is no A(t) in 3ord because sigma2_u in 3ord has the same role as 
  #   sigma2_a in exact method and for coding convenience.
  # lambda is NA if we want it to be estimated, as by default.
  #   If we want to fix it, we pass the value in the call of this function.
  sigma2_u <- ifelse(method.ssm == "3ord", 0, NA)
  sigma2_a <- NA
  sigma2_e <- NA
  
  
  # Z matrix
  Zt <- build_Zt(Z = Z,
                 S = S)
  
  # T matrix
  Tt <- build_Tt(method.ssm = method.ssm,
                 delta = delta,
                 lambda = lambda,
                 K = K)
  
  # R matrix
  Rt <- build_Rt(method.ssm = method.ssm,
                 K = K)
  
  # Q matrix
  Qt <- build_Qt(method.ssm = method.ssm,
                 sigma2_u = sigma2_u,
                 sigma2_a = sigma2_a,
                 lambda = lambda,
                 delta = delta,
                 kernel = kernel,
                 S = S)
  
  # H matrix
  Ht <- build_Ht(sigma2_e = sigma2_e,
                 Z = Z)
  
  # Initial values
  u0 <- build_u0(rates = rates,
                 S = S)
  a1 <- u0$a1
  P1 <- u0$P1
  
  # Latent states names.
  # In 3ord it would be more correct to call the 3rd parameter d2U, but
  #   I prefer to keep the same notation for coding convenience.
  state_names <- c(paste(c('U', 'dU', 'a'),
                         rep(0:K, each=3), sep=''))
  
  # Function to check well-definition of the model (for KFAS)
  checkfn <- function(model_gauss){
    # test positive semidefiniteness of H and Q
    !inherits(try(ldl(model_gauss$Q[,,1]), TRUE), 'try-error') &&
      !inherits(try(ldl(model_gauss$H[,,1]),TRUE),'try-error')
  }
  
  # Function to update the model components 
  # for a new set of hyperparameters (to run optimization for KFAS)
  updatefn <- function(pars, 
                       model_gauss,
                       method.ssm,
                       lambda,
                       S, 
                       Z, 
                       kernel, 
                       ages_max){
    
    K <- ncol(S)
    
    # Update the model with new hyperparameters
    lambda <- ifelse(is.na(lambda),
                     exp(pars['lambda']),
                     lambda)
    sigma2_u <- ifelse(method.ssm != "3ord",
                       exp(pars['sigma2_u']),
                       0)
    sigma2_a <- exp(pars['sigma2_a'])
    sigma2_e <- exp(pars['sigma2_e'])
    
    
    Qt <- build_Qt(method.ssm = method.ssm,
                   sigma2_u = sigma2_u,
                   sigma2_a = sigma2_a,
                   lambda = lambda,
                   delta = delta,
                   kernel = kernel,
                   S = S)
    model_gauss["Q", etas='custom'] <- Qt
    
    Tt <- build_Tt(method.ssm = method.ssm,
                   delta = delta,
                   lambda = lambda,
                   K = K)
    model_gauss["T"] <- Tt
    
    Ht <- build_Ht(sigma2_e = sigma2_e,
                   Z = Z)
    model_gauss["H"] <- Ht
    
    return(model_gauss)
  }
  
  # Model definition in KFAS
  model_gauss <- SSModel(log(rates) ~ -1 + 
                           SSMcustom(Z = Zt, T = Tt,
                                     R = Rt, Q = Qt,
                                     state_names = state_names,
                                     a1 = a1,
                                     P1 = P1,
                                     n = nrow(rates)), 
                         distribution = 'gaussian',
                         H = Ht)
  
  pars_tmp <- c('lambda', 'sigma2_u', 'sigma2_a', 'sigma2_e')
  pars <- pars_tmp[is.na(c(lambda, sigma2_u, sigma2_a, sigma2_e))]
  
  # Saving additional info 
  info <- append(info, list(S = S, 
                            K = K,
                            Z = Z,
                            ages_max = ages_max,
                            updatefn = updatefn,
                            checkfn = checkfn,
                            lambda = lambda,
                            pars = pars,
                            n_pars = length(pars),
                            lik = 'gaussian',
                            call = list('method.ssm' = method.ssm,
                                        'lambda.fixed' = !is.na(lambda))
                            )
                 )
  
  return(list('model' = model_gauss,
              'info' = info))
}



######################################################
# BSP model fit via MLE #############################
######################################################
#
# INPUT:
# - model: model returned from bsp.model
# - rep: how many time to repeat the optimization with
#        varying starting point
# - parallel: either to parallelize the "rep" runs
# - opt.method: optimization algorithm to use
# - maxcl: maximum number of cluster to use 
#          (if parallel = TRUE)
#
# RETURN:
# - List with fitted model and additional info
#

bsp.fit <- function(model, 
                    rep = 1, 
                    parallel = FALSE, 
                    opt.method = 'L-BFGS-B', 
                    maxcl = 10){
  
  info <- model$info
  call <- info$call
  
  cat("Fitting model",  call$method.ssm, 
      " with lambda", ifelse(call$lambda.fixed, "fixed", "free"), "\n")
  
  ## Function to minimise
  # logLik + inv-gamma prior on sigma2_u and sigma2_a
  fn_optim <- function(pars, 
                       model_gauss,
                       method.ssm,
                       lambda,
                       updatefn,
                       checkfn,
                       S, 
                       Z, 
                       kernel, 
                       ages_max) {
    
    
    # Update model with new pars value
    model_up <- updatefn(pars, 
                         model_gauss,
                         method.ssm,
                         lambda,
                         S, 
                         Z, 
                         kernel, 
                         ages_max)
    
    if(!checkfn(model_up)) stop('checkfn failed: model_up not valid. -> inside fn_optim')
    
    target_max <- logLik(model_up) + 
      dinvgamma(x = exp(pars['sigma2_a']), shape = 0.01, scale = 0.01, log = TRUE)
    
    # If not 3ord, add inv-gamma prior on sigma2_u. If 3ord then sigma2_u = 0,
    #   so no prior needed.
    if (method.ssm != "3ord"){
      target_max <- target_max + 
        dinvgamma(x = exp(pars['sigma2_u']), shape = 0.01, scale = 0.01, log = TRUE)
    }
      
    return(- target_max) # Because optim does minimization
  }
  
  # Starting points for the optimization algorithm
  starting_values <- lapply(1:rep, 
                            function(x) log(runif(info$n_pars, 
                                                  min = 1e-1, 
                                                  max = 1)))
  for (sv in seq_along(starting_values)){
    names(starting_values[[sv]]) <- info$pars
  }
  
  # Optimization
  if(parallel){
    registerDoParallel(cores = min(rep,
                                   maxcl))
    cl_fit <- makeCluster(min(rep,
                              maxcl),
                          type = ifelse(.Platform$OS.type == "unix",
                                        "FORK",
                                        "PSOCK"))
    
    clusterCall(cl_fit, function() {require(here); here('source', 'BSP.R')})
    
    fit_list <- parLapply(cl_fit,
                          starting_values,
                          function(init){
                            try(optim(par = init,
                                      fn = fn_optim,
                                      method = opt.method,
                                      control = list(maxit = 1e6),
                                      model_gauss = model$model,
                                      method.ssm = call$method.ssm,
                                      lambda = info$lambda,
                                      updatefn = info$updatefn,
                                      checkfn = info$checkfn,
                                      S = info$S,
                                      Z = info$Z,
                                      kernel = info$kernel,
                                      ages_max = info$ages_max), TRUE)
                          })
    stopCluster(cl_fit)
  } else{
    fit_list <- lapply(starting_values, 
                       function(init){
                         try(optim(par = init,
                                   fn = fn_optim,
                                   method = opt.method,
                                   control = list(maxit = 1e6),
                                   model_gauss = model$model,
                                   method.ssm = call$method.ssm,
                                   lambda = info$lambda,
                                   updatefn = info$updatefn,
                                   checkfn = info$checkfn,
                                   S = info$S,
                                   Z = info$Z,
                                   kernel = info$kernel,
                                   ages_max = info$ages_max), FALSE)
                       })
  }
  
  # Throw away possible initializations that led to errors
  fit_list_clean <- discard(fit_list, . %>% 
                              inherits(., 'try-error'))
  print(paste('Number of failed optim attempt:', 
              length(fit_list) - length(fit_list_clean)))
  if(length(fit_list_clean) == 0){
    print('All optimization attempts failed..')
    return(-1)
  }
  
  # Keep track of the optimization history -> internal use
  history <- as_tibble(sapply(fit_list_clean, . %>% `$`(., par))) %>%
    bind_rows(as_tibble(matrix(sapply(fit_list_clean, . %>% `$`(., value)), nrow = 1))) %>%
    mutate(names = c(info$pars, 'neg_p_loglik'))
  
  if (info$call$lambda.fixed){ 
    history <- history %>%
      bind_rows(matrix(c(rep(info$lambda, rep)), 
                       nrow = 1) %>% 
                  as.data.frame() %>%
                  mutate(names = "lambda (fixed)"))
  }
  
  if (info$call$method.ssm == '3ord'){
    history <- history %>%
      bind_rows(matrix(c(rep(NA, rep)), 
                       nrow = 1) %>% 
                  as.data.frame() %>%
                  mutate(names = "sigma_u (not estimated)"))
  }
  
  # Select best parameters found
  best_fit <- which.min(sapply(fit_list_clean, . %>% 
                                 `$`(.,value)))
  fit_bin <- fit_list_clean[[best_fit]]
  
  info <- append(info, list(optim = fit_bin,
                            rep = rep,
                            call = 'bsp.fit',
                            opt.method = opt.method,
                            failed_optim = length(fit_list) - length(fit_list_clean)))
  
  
  return(list('fit' = info$updatefn(pars = fit_bin$par,
                                    model_gauss = model$model, 
                                    method.ssm = call$method.ssm,
                                    lambda = info$lambda,
                                    S = info$S,
                                    Z = info$Z,
                                    kernel = info$kernel,
                                    ages_max = info$ages_max),
              'info' = info,
              'history' = history))
}
