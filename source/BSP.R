require(magrittr)
require(KFAS)
require(splines2)
require(invgamma)
require(doParallel)

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

bsp.model <- function(rates, 
                      delta, 
                      age_knots,
                      kernel) {
  info <- list()
  info <- append(info, list(data.nrow = nrow(rates),
                            data.ncol = ncol(rates),
                            age_knots = age_knots,
                            kernel = kernel,
                            delta = delta))

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
  
  #### State-space model
  #
  # u_t = (u_t1, ..., u_tK) stays for 
  # beta_t = (beta_t1, ..., beta_tK) in the paper.
  #
  # Model structure in KFAS:
  # log_m_t = Zt*beta_t + N(0,Ht)
  # beta_t+1 = Tt*beta_t + Rt*N(0,Qt)
  
  # NA are hyperparameters to estimate via MLE
  lambda <- NA
  sigma2_u <- NA
  sigma2_a <- NA
  sigma2_e <- NA

  # Z matrix
  Zt <- matrix(0, nrow = Z, ncol = 3*(K+1))
  Zt[-1, seq(4, 3*(K+1), by=3)] <- S
  Zt[1,1] <- 1
  
  # T matrix
  T_tmp <- diag(1,3)
  T_tmp[1,2] <- T_tmp[2,3] <- delta*lambda
  T_tmp[1,3] <- delta^2/2*lambda^2
  Tt <- kronecker(diag(1,K+1), T_tmp)
  
  # R matrix
  Rt <- diag(1,3*(K+1))
  
  # Creating covariance matrix Q
  ages_max <- c(0, apply(S, 2, which.max))
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
  
  # Observational variance
  sigma2_e <- NA
  Ht <- diag(sigma2_e, Z)
  
  # Initial values u_0 = a1 set at frequentist estimate 
  Eta <- log(rates[,-1])
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
  
  # Latent states name
  # In the paper: U -> beta, dU -> dBeta, a -> a
  state_names <- c(paste(c('U', 'dU', 'a'),
                         rep(0:K, each=3), sep=''))
  
  # Function to check well-definition of the model (for KFAS)
  checkfn <- function(model){
    # test positive semidefiniteness of Q and positivity of H
    (model["Q"] > 0) &&
      !inherits(try(ldl(model$Q[,,1]), TRUE), 'try-error')
    (model["H"] > 0) &&
      !inherits(try(ldl(model$H[,,1]),TRUE),'try-error')
  }
  
  # Function to update the model components 
  # for a new set of hyperparameters (to run optimization for KFAS)
  updatefn <- function(pars, 
                       model, 
                       K, 
                       Z, 
                       kernel, 
                       ages_max){
    lambda <- exp(pars[1])
    sigma2_u <- exp(pars[2])
    sigma2_a <- exp(pars[3])
    sigma2_e <- exp(pars[4])
    Q_tmp <- matrix(NA, 3, 3)
    Q_tmp[1,1] <- sigma2_u/3*delta^3*lambda^2 + sigma2_a/20*delta^5*lambda^4
    Q_tmp[1,2] <- Q_tmp[2,1] <- sigma2_u/2*delta^2*lambda + sigma2_a/8*delta^4*lambda^3
    Q_tmp[1,3] <- Q_tmp[3,1] <- sigma2_a/6*delta^3*lambda^2
    Q_tmp[2,2] <- sigma2_u*delta + sigma2_a/3*delta^3*lambda^2
    Q_tmp[2,3] <- Q_tmp[3,2] <- sigma2_a/2*delta^2*lambda
    Q_tmp[3,3] <- sigma2_a*delta
    Qt <- kronecker(diag(1,K+1), Q_tmp)
    for(I in 1:(K+1-1)){
      for(J in (I+1):(K+1)){
        rho_u <- kernel(x = abs(ages_max[J]-ages_max[I])) 
        rho_a <- 0
        Q_tmp <- matrix(NA, 3, 3)
        Q_tmp[1,1] <- sigma2_u/3*rho_u*delta^3*lambda^2 + sigma2_a/20*rho_a*delta^5*lambda^4
        Q_tmp[1,2] <- Q_tmp[2,1] <- sigma2_u/2*rho_u*delta^2*lambda + 
          sigma2_a/8*rho_a*delta^4*lambda^3
        Q_tmp[1,3] <- Q_tmp[3,1] <- sigma2_a/6*rho_a*delta^3*lambda^2
        Q_tmp[2,2] <- sigma2_u*rho_u*delta + sigma2_a/3*rho_a*delta^3*lambda^2
        Q_tmp[2,3] <- Q_tmp[3,2] <- sigma2_a/2*rho_a*delta^2*lambda
        Q_tmp[3,3] <- sigma2_a*rho_a*delta
        Qt[(1+(I-1)*3):(I*3),(1+(J-1)*3):(J*3)] <- 
          Qt[(1+(J-1)*3):(J*3),(1+(I-1)*3):(I*3)] <- Q_tmp
      }
    }
    model["Q", etas='custom'] <- Qt
    T_tmp <- diag(1,3)
    T_tmp[1,2] <- T_tmp[2,3] <- delta*lambda
    T_tmp[1,3] <- delta^2/2*lambda^2
    Tt <- kronecker(diag(1,K+1), T_tmp)
    model["T"] <- Tt
    model["H"] <- diag(exp(pars[4]), Z)
    model
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
  
  # Saving additional info 
  info <- append(info, list(K = K,
                            Z = Z,
                            ages_max = ages_max,
                            updatefn = updatefn,
                            checkfn = checkfn,
                            n_pars = 4, # lambda, sigma2_e, sigma2_u, sigma2_e
                            lik = 'gaussian',
                            call = 'bsp.model'))
  
  return(list('model' = model_gauss,
              'info' = info))
}

ngp.model <- function(rates, 
                      delta) {
  info <- list()
  info <- append(info, list(data.nrow = nrow(rates),
                            data.ncol = ncol(rates),
                            delta = delta))
  
  # Maximum age
  Z <- ncol(rates)
  ages <- as.numeric(colnames(rates))
  ages_no0 <- ages[-1]
  
  # One nGP for each age
  K <- length(ages_no0) # 100
  
  #### State-space model
  #
  # u_t = (u_t1, ..., u_tK) stays for 
  # beta_t = (beta_t1, ..., beta_tK) in the paper.
  #
  # Model structure in KFAS:
  # log_m_t = Zt*beta_t + N(0,Ht)
  # beta_t+1 = Tt*beta_t + Rt*N(0,Qt)
  
  # NA are hyperparameters to estimate via MLE
  lambda <- 1
  sigma2_u <- NA
  sigma2_a <- NA
  sigma2_e <- NA
  
  # Z matrix
  Zt <- matrix(0, nrow = Z, ncol = 3*(K+1))
  Zt[, seq(1, 3*(K+1), by=3)] <- diag(1, nrow = Z, ncol = (K+1))
  
  # T matrix
  T_tmp <- diag(1,3)
  T_tmp[1,2] <- T_tmp[2,3] <- delta*lambda
  T_tmp[1,3] <- delta^2/2*lambda^2
  Tt <- kronecker(diag(1,K+1), T_tmp)
  
  # R matrix
  Rt <- diag(1,3*(K+1))
  
  # Creating covariance matrix Q
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
  # for(I in 1:(K+1-1)){
  #   for(J in (I+1):(K+1)){
  #     rho_u <- 0 # No correlation structure
  #     rho_a <- 0 
  #     Q_tmp[1,1] <- sigma2_u/3 * rho_u * delta^3 * lambda^2 + 
  #       sigma2_a/20 * rho_a * delta^5 * lambda^4
  #     Q_tmp[1,2] <- 
  #       Q_tmp[2,1] <- sigma2_u/2 * rho_u*delta^2 * lambda + 
  #       sigma2_a/8 * rho_a*delta^4 * lambda^3
  #     Q_tmp[1,3] <- 
  #       Q_tmp[3,1] <- sigma2_a/6 * rho_a * delta^3 * lambda^2
  #     Q_tmp[2,2] <- sigma2_u * rho_u * delta + 
  #       sigma2_a/3 * rho_a*delta^3 * lambda^2
  #     Q_tmp[2,3] <- 
  #       Q_tmp[3,2] <- sigma2_a/2 * rho_a * delta^2 * lambda
  #     Q_tmp[3,3] <- sigma2_a*rho_a*delta
  #     Qt[(1+(I-1)*3):(I*3),(1+(J-1)*3):(J*3)] <- 
  #       Qt[(1+(J-1)*3):(J*3),(1+(I-1)*3):(I*3)] <- Q_tmp
  #   }
  # }
  
  # Observational variance
  sigma2_e <- NA
  Ht <- diag(sigma2_e, Z)
  
  # Initial values u_0 = a1 set at observed rates
  tmp_a1 <- rep(0, 3*(K+1))
  tmp_a1[seq(1, 3*(K+1), by=3)] <- log(rates[1,])
  a1 <- matrix(tmp_a1, 3*(K+1), 1)
  
  #Initial variance var(u_0) = P1
  P1 <- diag(10, 3*(K+1))
  
  # Latent states name
  # In the paper: U -> beta, dU -> dBeta, a -> a
  state_names <- c(paste(c('U', 'dU', 'a'),
                         rep(0:K, each=3), sep=''))
  
  # Function to check well-definition of the model (for KFAS)
  checkfn <- function(model){
    # test positive semidefiniteness of Q and positivity of H
    (model["Q"] > 0) &&
      !inherits(try(ldl(model$Q[,,1]), TRUE), 'try-error')
    (model["H"] > 0) &&
      !inherits(try(ldl(model$H[,,1]),TRUE),'try-error')
  }
  
  # Function to update the model components 
  # for a new set of hyperparameters (to run optimization for KFAS)
  updatefn <- function(pars, 
                       model, 
                       K, 
                       Z){
    lambda <- 1
    sigma2_u <- exp(pars[1])
    sigma2_a <- exp(pars[2])
    sigma2_e <- exp(pars[3])
    Q_tmp <- matrix(NA, 3, 3)
    Q_tmp[1,1] <- sigma2_u/3*delta^3*lambda^2 + sigma2_a/20*delta^5*lambda^4
    Q_tmp[1,2] <- Q_tmp[2,1] <- sigma2_u/2*delta^2*lambda + sigma2_a/8*delta^4*lambda^3
    Q_tmp[1,3] <- Q_tmp[3,1] <- sigma2_a/6*delta^3*lambda^2
    Q_tmp[2,2] <- sigma2_u*delta + sigma2_a/3*delta^3*lambda^2
    Q_tmp[2,3] <- Q_tmp[3,2] <- sigma2_a/2*delta^2*lambda
    Q_tmp[3,3] <- sigma2_a*delta
    Qt <- kronecker(diag(1,K+1), Q_tmp)
    # for(I in 1:(K+1-1)){
    #   for(J in (I+1):(K+1)){
    #     rho_u <- 0
    #     rho_a <- 0
    #     Q_tmp <- matrix(NA, 3, 3)
    #     Q_tmp[1,1] <- sigma2_u/3*rho_u*delta^3*lambda^2 + sigma2_a/20*rho_a*delta^5*lambda^4
    #     Q_tmp[1,2] <- Q_tmp[2,1] <- sigma2_u/2*rho_u*delta^2*lambda + 
    #       sigma2_a/8*rho_a*delta^4*lambda^3
    #     Q_tmp[1,3] <- Q_tmp[3,1] <- sigma2_a/6*rho_a*delta^3*lambda^2
    #     Q_tmp[2,2] <- sigma2_u*rho_u*delta + sigma2_a/3*rho_a*delta^3*lambda^2
    #     Q_tmp[2,3] <- Q_tmp[3,2] <- sigma2_a/2*rho_a*delta^2*lambda
    #     Q_tmp[3,3] <- sigma2_a*rho_a*delta
    #     Qt[(1+(I-1)*3):(I*3),(1+(J-1)*3):(J*3)] <- 
    #       Qt[(1+(J-1)*3):(J*3),(1+(I-1)*3):(I*3)] <- Q_tmp
    #   }
    # }
    model["Q", etas='custom'] <- Qt
    T_tmp <- diag(1,3)
    T_tmp[1,2] <- T_tmp[2,3] <- delta*lambda
    T_tmp[1,3] <- delta^2/2*lambda^2
    Tt <- kronecker(diag(1,K+1), T_tmp)
    model["T"] <- Tt
    model["H"] <- diag(sigma2_e, Z)
    model
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
  
  # Saving additional info 
  info <- append(info, list(K = K,
                            Z = Z,
                            updatefn = updatefn,
                            checkfn = checkfn,
                            n_pars = 3, # sigma2_e, sigma2_u, sigma2_e
                            lik = 'gaussian',
                            call = 'ngp.model'))
  
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
# - method: optimization algorithm to use
# - maxcl: maximum number of cluster to use 
#          (if parallel = TRUE)
#
# RETURN:
# - List with fitted model and additional info
#

bsp.fit <- function(model, 
                     rep = 1, 
                     parallel = FALSE, 
                     method = 'L-BFGS-B', 
                     maxcl = 30){
  model_gauss <- model$model
  info <- model$info
  updatefn <- info$updatefn
  checkfn <- info$checkfn
  
  ## Function to minimise
  # logLik + inv-gamma prior on sigma2_u and sigma2_a
  fn_optim <- function(pars, 
                       model, 
                       updatefn, 
                       checkfn, 
                       K, 
                       Z, 
                       kernel, 
                       ages_max) {
    # Update model with new pars value
    model_up <- updatefn(pars, 
                         model, 
                         K, 
                         Z, 
                         kernel, 
                         ages_max)
    
    if(!checkfn(model_up)) stop('checkfn failed: model_up not valid. -> inside fn_optim')
    
    target_max <- logLik(model_up) + 
      dinvgamma(x = exp(pars[2]), shape = 0.01, scale = 0.01, log = TRUE) +
      dinvgamma(x = exp(pars[3]), shape = 0.01, scale = 0.01, log = TRUE) 
    return(- target_max) # Because optim does minimization
  }
  
  # Starting points for the optimization algorithm
  starting_values <- lapply(1:rep, 
                            function(x) log(runif(info$n_pars, 
                                                  min = 1e-1, 
                                                  max = 1)))
  # Optimization
  if(parallel){
    registerDoParallel(cores = min(rep, 
                                   maxcl))
    cl_fit <- makeCluster(min(rep, 
                              maxcl), 
                          type = "FORK")
    fit_list <- parLapply(cl_fit,
                          starting_values, 
                          function(init){
                            try(optim(par = init,
                                      fn = fn_optim,
                                      method = method,
                                      control = list(maxit = 1e6),
                                      model = model_gauss,
                                      updatefn = updatefn,
                                      checkfn = checkfn,
                                      K = info$K,
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
                                   method = method,
                                   control = list(maxit = 1e6),
                                   model = model_gauss,
                                   updatefn = updatefn,
                                   checkfn = checkfn,
                                   K = info$K,
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
    mutate(names = c('lambda', 'sigma2_u', 'sigma2_a', 'sigma_e', 'neg_p_loglik'))
  
  # Select best parameters found
  best_fit <- which.min(sapply(fit_list_clean, . %>% 
                                 `$`(.,value)))
  fit_bin <- fit_list_clean[[best_fit]]
  
  info <- append(info, list(optim = fit_bin,
                            rep = rep,
                            call = 'bsp.fit',
                            method = method,
                            failed_optim = length(fit_list) - length(fit_list_clean)))
  
  return(list('fit' = updatefn(fit_bin$par,
                               model_gauss, 
                               info$K, 
                               info$Z, 
                               info$kernel, 
                               info$ages_max),
              'info' = info,
              'history' = history))
}

ngp.fit <- function(model, 
                    rep = 1, 
                    parallel = FALSE, 
                    method = 'L-BFGS-B', 
                    maxcl = 30){
  model_gauss <- model$model
  info <- model$info
  updatefn <- info$updatefn
  checkfn <- info$checkfn
  
  ## Function to minimise
  # logLik + inv-gamma prior on sigma2_u and sigma2_a
  fn_optim <- function(pars, 
                       model, 
                       updatefn, 
                       checkfn, 
                       K, 
                       Z) {
    # Update model with new pars value
    model_up <- updatefn(pars, 
                         model, 
                         K, 
                         Z)
    
    if(!checkfn(model_up)) stop('checkfn failed: model_up not valid. -> inside fn_optim')
    
    target_max <- logLik(model_up) + 
      dinvgamma(x = exp(pars[2]), shape = 0.01, scale = 0.01, log = TRUE) +
      dinvgamma(x = exp(pars[3]), shape = 0.01, scale = 0.01, log = TRUE) 
    return(- target_max) # Because optim does minimization
  }
  
  # Starting points for the optimization algorithm
  starting_values <- lapply(1:rep, 
                            function(x) log(runif(info$n_pars, 
                                                  min = 1e-1, 
                                                  max = 1)))
  # Optimization
  if(parallel){
    registerDoParallel(cores = min(rep, 
                                   maxcl))
    cl_fit <- makeCluster(min(rep, 
                              maxcl), 
                          type = "FORK")
    fit_list <- parLapply(cl_fit,
                          starting_values, 
                          function(init){
                            try(optim(par = init,
                                      fn = fn_optim,
                                      method = method,
                                      control = list(maxit = 1e6),
                                      model = model_gauss,
                                      updatefn = updatefn,
                                      checkfn = checkfn,
                                      K = info$K,
                                      Z = info$Z), TRUE)
                          })
    stopCluster(cl_fit)
  } else{
    fit_list <- lapply(starting_values, 
                       function(init){
                         try(optim(par = init,
                                   fn = fn_optim,
                                   method = method,
                                   control = list(maxit = 1e6),
                                   model = model_gauss,
                                   updatefn = updatefn,
                                   checkfn = checkfn,
                                   K = info$K,
                                   Z = info$Z), FALSE)
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
    mutate(names = c('sigma2_u', 'sigma2_a', 'sigma_e', 'neg_p_loglik'))
  
  # Select best parameters found
  best_fit <- which.min(sapply(fit_list_clean, . %>% 
                                 `$`(.,value)))
  fit_bin <- fit_list_clean[[best_fit]]
  
  info <- append(info, list(optim = fit_bin,
                            rep = rep,
                            call = 'ngp.fit',
                            method = method,
                            failed_optim = length(fit_list) - length(fit_list_clean)))
  
  return(list('fit' = updatefn(fit_bin$par,
                               model_gauss, 
                               info$K, 
                               info$Z),
              'info' = info,
              'history' = history))
}

