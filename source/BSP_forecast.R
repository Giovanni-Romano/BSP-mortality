# This script is mostly designed for the increasing sample size
# forecasting esperiment implemented in the paper.
# However, user can find useful functions optim_rwd()
# and fitandforecast() as they can be used 
# in more general settings.

require(mvtnorm)

###############################################
## Fit and h-step-ahead forecast  #############
###############################################
#
# INPUT:
#   - model: BSP model to estimate (result of bsp.model)
#   - h: how many step-ahead to forecast
#   - Z: maximum age
#   - rep: how many time to repeat optimization
#   - n_for: how many years to use in the trainig data
#           to estimate the rw+drift model
#   - method: optimization method
#   - parallel: parallelization
#   - maxcl: maximum number of cluster ofr the parallelization
#
# OUTPUT:
#   - pred: forecasts
#

fitandforecast_uq <- function(model, 
                               h, 
                               Z, 
                               rep = 5, 
                               n_for, 
                               opt.method = 'Nelder-Mead', 
                               parallel = FALSE, 
                               maxcl = 10) 
{
  fit <- bsp.fit(model, 
                 rep = rep, 
                 opt.method = opt.method, 
                 parallel = parallel, 
                 maxcl = maxcl)
  K <- fit$info$K
  Z <- model$info$Z
  ## Fitting random-walk with drift for prediction
  # Simulation from the smoothing distribution
  # Smoothing distribution
  smooth <- KFS(fit$fit, smoothing = c('state','signal'))
  sim_states <- simulateSSM(fit$fit, 
                            type = 'states',
                            filtered = FALSE,
                            conditional = TRUE,
                            antithetics = TRUE,
                            nsim = 100)
  sim_last25 <- tail(sim_states, 
                     n = c(n_for,3*(K+1),dim(sim_states)[3]))
  sim_last25_B25 <- tail(sim_states, 
                         n = c(2*n_for,3*(K+1),dim(sim_states)[3]))[1:n_for,,]
  U <- apply(sim_last25[,paste("U",0:19,sep=""),],
             MARGIN = c(1,2), mean)

  U_m_B25 <- apply(sim_states[train - 2*n_for, paste('U', 0:K, sep = ''),],
                   MARGIN = 1, mean)
  dU <- apply(sim_last25[,paste("dU",0:19,sep=""),],
              MARGIN = c(1,2), mean) %>%
    apply(MARGIN = 2, median)
  dU_B25 <- apply(sim_last25_B25[,paste("dU",0:19,sep=""),],
                  MARGIN = c(1,2), mean) %>%
    apply(MARGIN = 2, median)
  drift_median_var <- apply(sim_last25[,paste("dU",0:19,sep=""),],
                            MARGIN = c(2,3), median) %>%
    apply(MARGIN = 1, var)
  drift_median_var_B25 <- apply(sim_last25_B25[,paste("dU",0:19,sep=""),],
                                MARGIN = c(2,3), median) %>%
    apply(MARGIN = 1, var)

  ## Creating random walk + dirft model
  Tt_pred <- fit$fit$T[1,2,1] # lambda*delta
  kernel <- fit$info$kernel
  ages_max <- fit$info$ages_max
  
  ## Composing rw + drift with estimated sigma2_u
  sigma2_u <- NA
  Qt_pred <- diag(2*(K+1))
  diag(Qt_pred)[1:(K+1)] <- sigma2_u
  diag(Qt_pred)[-c(1:(K+1))] <- NA 
  for(I in 1:(K+1-1)){
    for(J in (I+1):(K+1)){
      rho_u <- kernel(x = abs(ages_max[J]-ages_max[I]))
      Qt_pred[I,J] <- Qt_pred[J,I] <- sigma2_u * rho_u
    }
  }
  Rt_pred <- diag(1, 2*(K+1), 2*(K+1))
  Tt_pred_mat <- diag(1, 2*(K+1), 2*(K+1))
  Tt_pred_mat[1:(K+1),(K+2):(2*(K+1))] <- diag(1, K+1, K+1)
  
  Zt_pred_ <- fit$fit$Z[, seq(1, 3*(K+1), by = 3), ]
  Zt_pred <- matrix(0, Z, 2*(K+1))
  Zt_pred[, 1:(K+1)] <- Zt_pred_
  a1_pred <- matrix(c(U_m_B25 + Tt_pred*dU_B25,
                      Tt_pred*dU_B25), 
                    2*(K+1), 1) # last year of training

  P1_pred <- matrix(0, 2*(K+1), 2*(K+1))
  P1_pred[1:(K+1), 1:(K+1)] <- smooth$P[seq(1, 3*(K+1), by = 3),seq(1, 3*(K+1), by = 3),train-n_for+1] # +
   # diag(Tt_pred^2 * drift_median_var_B25)
  diag(P1_pred[-c(1:(K+1)), -c(1:(K+1))]) <- Tt_pred^2 * drift_median_var_B25

  Ht_pred <- diag(NA, Z)
  state_names <- c(paste('U', 0:K, sep = ''),
                   paste('drift', 0:K, sep = ''))
  
  updatefn <- function(pars, 
                       model,
                       Z,
                       delta_lambda_sq){
    sigma2_u <- exp(pars[1])
    drift_var <- exp(pars[2])
    Qt_pred <- diag(2*(K+1))
    diag(Qt_pred)[1:(K+1)] <- sigma2_u
    diag(Qt_pred)[-c(1:(K+1))] <- delta_lambda_sq*drift_var
    for(I in 1:(K+1-1)){
      for(J in (I+1):(K+1)){
        rho_u <- kernel(x = abs(ages_max[J]-ages_max[I]))
        Qt_pred[I,J] <- Qt_pred[J,I] <- sigma2_u * rho_u
      }
    }
    model["Q"] <- Qt_pred
    model["H"] <- diag(exp(pars[3]), Z)
    model
  }

  rwd_gauss <- SSModel(tail(model$model$y, n = n_for, k = 1) ~ -1 + 
                         SSMcustom(Z = Zt_pred, T = Tt_pred_mat,
                                   R = Rt_pred, Q = Qt_pred,
                                   state_names = state_names,
                                   a1 = a1_pred,
                                   P1 = P1_pred,
                                   # P1inf = diag(1, 2*(K+1), 2*(K+1)),
                                   n = nrow(tail(model$model$y, n = n_for, k = 1))), 
                       distribution = 'gaussian',
                       H = Ht_pred)
  starting_values <- lapply(1:3, 
                            function(x) log(runif(3, min = 1e-4, max = 2)))
  rwd_kfas_fit_list <- lapply(starting_values,
                              function(init){
                                fitSSM(rwd_gauss,
                                       inits = init,
                                       method = opt.method,
                                       updatefn = updatefn,
                                       checkfn = model$info$checkfn,
                                       update_args = list(Z = Z,
                                                          delta_lambda_sq = Tt_pred^2))
                              })
  best_fit <- which.min(sapply(rwd_kfas_fit_list, . %>% 
                                 `$`(.,optim.out) %>%
                                 `$`(.,value)))
  rwd_kfas_fit <- rwd_kfas_fit_list[[best_fit]]

  rwd_gauss_fit <- rwd_kfas_fit$model
  rwd_smooth <- KFS(rwd_gauss_fit,
                    smoothing = c('state'))
  V_rwd <- rwd_smooth$V[,,n_for]
  Qt_pred <- rwd_gauss_fit$Q[,,1] 
  H_rwd <- rwd_gauss_fit$H[,,1] 
  # Prediction h-step ahead
  a1_pred <- tail(U, n = 1 ,k = 1) # Now a1 is last observed year
  U_pred <- matrix(NA, nrow = h, ncol = K+1)
  varU_ext_pred <- array(NA, dim = c(2*(K+1), 2*(K+1), h))
  f_pred <- matrix(NA, nrow = h, ncol = Z)
  varf_pred <- array(NA, dim = c(Z, Z, h))
  for(tt in 1:h){
    if(tt == 1){
      U_pred[1,] <- as.numeric(a1_pred[1:(K+1)]) + Tt_pred%*%dU
      varU_ext_pred[,,1] <- Tt_pred_mat %*% V_rwd %*% t(Tt_pred_mat) + Qt_pred
      f_pred[1,] <- Zt_pred_%*%U_pred[1,]
      varf_pred[,,1] <- H_rwd + 
        Zt_pred %*% varU_ext_pred[,,1] %*% t(Zt_pred)
    } else {
      U_pred[tt,] <- U_pred[tt-1,] + Tt_pred%*%dU
      varU_ext_pred[,,tt] <- Tt_pred_mat %*% varU_ext_pred[,,tt-1] %*% Tt_pred_mat + Qt_pred 
      f_pred[tt,] <- Zt_pred_%*%U_pred[tt,]
      varf_pred[,,tt] <- H_rwd + 
        Zt_pred %*% varU_ext_pred[,,tt] %*% t(Zt_pred)
    }
  }
  pred <- list()
  for(zz in 1:Z){
    pred <- c(pred,
              list(tibble(fit = f_pred[,zz],
                          upr = qnorm(p = 0.975,
                                      mean = f_pred[,zz],
                                      sd = sqrt(varf_pred[zz,zz,])),
                          lwr = qnorm(p = 0.025,
                                      mean = f_pred[,zz],
                                      sd = sqrt(varf_pred[zz,zz,])))))
  }
  names(pred) <- as.character(0:(Z-1))
  return(pred)
}



###############################################
## Helper function to fit with increasing     #
## sample size (interal use)                  #
###############################################
rolling_uq <- function(cg,
                       method.ssm,
                       lambda = NA,
                       n_for = 25, 
                       parallel = FALSE,
                       maxcl = 10)
{
  # Loading the data corrisponding to country cg
  country <- sub("_.*", "", cg)
  gender <- sub(".*_", "", cg)
  print(paste('Doing', country, gender))
  Y <- eval(parse(text = paste('Y', country, gender, sep = '_')))
  N <- eval(parse(text = paste('N', country, gender, sep = '_')))
  Tmax <- nrow(Y)
  Z <- ncol(Y)
  
  # Creating increasing sample size datasets
  data_list <- lapply(train:(Tmax - h_step),
                      function(t){
                        Y[1:t,]/N[1:t,]
                      })
  model_list <- lapply(data_list, . %>% bsp.model(method.ssm = method.ssm,
                                                  lambda = lambda,
                                                  rates = .,
                                                  delta = delta,
                                                  age_knots = age_knots,
                                                  kernel = matern_kernel))
  # Fit and forecast
  fitandfor_list <- lapply(model_list,
                           . %>% fitandforecast_uq(.,
                                                   h = h_step,
                                                   Z = Z,
                                                   rep = rep,
                                                   n_for = n_for,
                                                   opt.method = 'Nelder-Mead',
                                                   parallel = parallel,
                                                   maxcl = 30))
  # Postprocessing of forecast list
  forecast_list <- try(fitandfor_list %>%
                         modify(. %>% imap_dfr(.f = ~ (as_tibble(.) %>%
                                                         mutate(h_ahead = 1:h_step)),
                                               .id = 'age')) %>%
                         imap_dfr(.f = ~ ., .id='t') %>%
                         mutate_at(vars(t), as.numeric) %>%
                         mutate(country = country,
                                gender = gender))
  return(forecast_list)
}