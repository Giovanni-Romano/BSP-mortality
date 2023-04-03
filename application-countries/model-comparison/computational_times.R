require(here)
require(tidyverse)
require(KFAS)
require(StMoMo)
require(MortalitySmooth) 
require(svcm)
require(MortalityForecast)
require(bench)

rm(list=ls())

set.seed(132342)

output_collector <- list()
output_collector <- append(output_collector, list(info=Sys.time()))

source(here('source','BSP.R'))
source(here('source','setup.R'))
source(here('source','BSP_forecast.R'))
source(here('application-countries',
            'model-comparison',
            'CP',
            'SmoothConstrainedMortalityForecasting_Functions.R'))
source(here('application-countries',
            'model-comparison',
            'CP',
            'SmoothConstrainedMortalityForecasting_LifeTableFunctions.R'))

forecast_v2 <- function()
{
  K <- fit$info$K
  Z <- model$info$Z
  ## Fitting random-walk with drif for prediction
  # Smoothing distribution
  smooth <- KFS(fit$fit, smoothing = c('state','signal'))
  # Extracting info on U and dU from BSP model
  U_m <- smooth$alphahat[, paste('U', 0:K, sep = '')][train - n_for + 1,]
  U <- smooth$alphahat[, paste('U', 0:K, sep = '')] %>% 
    tail(n = n_for, k = 1) # selecting last n_for years
  dU <- apply(smooth$alphahat[, paste('dU', 0:K, sep = '')] %>%
                tail(n = n_for, k = 1), # selecting last n_for years
              MARGIN = 2, median)
  ## Creating random walk + dirft model
  Tt_pred <- fit$fit$T[1,2,1] # lambda*delta
  kernel <- fit$info$kernel
  ages_max <- fit$info$ages_max
  # Estimation of rw + drift
  rep_pred <- 1
  starting_values <- log(runif(rep_pred, min = 1e-2, max = 2))
  fit_rw_list <- lapply(starting_values, 
                        function(init){
                          try(optim(par = init,
                                    fn = optim_rwd,
                                    method = "L-BFGS-B",
                                    control = list(maxit = 1e6),
                                    u = U, 
                                    mu0 = U_m,
                                    Tt = Tt_pred,
                                    drift = dU,
                                    kernel = fit$info$kernel,
                                    ages_max = fit$info$ages_max,
                                    K = K), TRUE)
                        })
  # Checking failed optim
  fit_rw_list_clean <- discard(fit_rw_list, 
                               . %>% inherits(., 'try-error'))
  print(paste('Number of failed optim attempt:', 
              length(fit_rw_list) - length(fit_rw_list_clean)))
  if(length(fit_rw_list_clean) == 0){
    print('All optimization attempts failed..')
    return(-1)
  }
  # Extracting best fit
  best_fit_rw <- which.min(sapply(fit_rw_list_clean, . %>% `$`(.,value)))
  fit_rw <- fit_rw_list_clean[[best_fit_rw]]
  
  ## Composing rw + drift with estimated sigma2_u
  sigma2_u <- exp(fit_rw$par)
  Qt_pred <- diag(K+1) * sigma2_u
  for(I in 1:(K+1-1)){
    for(J in (I+1):(K+1)){
      rho_u <- kernel(x = abs(ages_max[J]-ages_max[I]))
      Qt_pred[I,J] <- Qt_pred[J,I] <- sigma2_u * rho_u
    }
  }
  Rt_pred <- matrix(0, 2*(K+1), K+1)
  Rt_pred[1:(K+1),] <- diag(1, K+1, K+1)
  Tt_pred_mat <- diag(1, 2*(K+1), 2*(K+1))
  Tt_pred_mat[1:(K+1),(K+2):(2*(K+1))] <- diag(1, K+1, K+1)
  
  Zt_pred_ <- fit$fit$Z[, seq(1, 3*(K+1), by = 3), ]
  Zt_pred <- matrix(0, Z, 2*(K+1))
  Zt_pred[, 1:(K+1)] <- Zt_pred_
  a1_pred <- matrix(c(tail(U, n = 1 ,k = 1),
                      Tt_pred*dU), 
                    2*(K+1), 1) # last year of training
  P1_pred_ <- smooth$V %>%
    tail(n = c(3*(K+1), 3*(K+1), 1)) %>% # selecting last year 
    # extracting only var(u) components
    `[`(seq(1, 3*(K+1), by = 3), seq(1, 3*(K+1), by = 3), ) 
  P1_pred <- diag(1, 2*(K+1), 2*(K+1))
  P1_pred[1:(K+1), 1:(K+1)] <- P1_pred_
  Ht_pred <- diag(NA, Z)
  state_names <- c(paste('U', 0:K, sep = ''),
                   paste('drift', 0:K, sep = ''))
  
  ## Fit rw-error sigma_e to data 
  # Use KFAS!
  updatefn <- function(pars, 
                       model,
                       Z){
    model["H"] <- diag(exp(pars[1]), Z)
    model
  }
  rwd_gauss <- SSModel(tail(data, n = n_for, k = 1) ~ -1 + 
                         SSMcustom(Z = Zt_pred, T = Tt_pred_mat,
                                   R = Rt_pred, Q = Qt_pred,
                                   state_names = state_names,
                                   a1 = a1_pred,
                                   P1 = P1_pred,
                                   n = nrow(tail(model$model$y, n = n_for, k = 1))), 
                       distribution = 'gaussian',
                       H = Ht_pred)
  rwd_kfas_fit <- fitSSM(rwd_gauss, 
                         inits = log(runif(1, min = 1e-2, max = 2)), 
                         method = "L-BFGS-B",
                         updatefn = updatefn,
                         checkfn = model$info$checkfn,
                         update_args = list(Z = Z))
  pred <- predict(rwd_kfas_fit$model,
                  interval = 'prediction',
                  newdata = SSModel(matrix(NA, nrow = h, ncol = Z) ~ -1 +
                                      SSMcustom(Z = rwd_kfas_fit$model$Z, T = rwd_kfas_fit$model$T,
                                                R = rwd_kfas_fit$model$R, Q = rwd_kfas_fit$model$Q,
                                                state_names = rownames(rwd_kfas_fit$model$a1),
                                                n = h), # total number of years, dimensionality check
                                    distribution = 'gaussian',
                                    H = rwd_kfas_fit$model$H),
                  level = 0.95,
                  type = 'response')
  return(pred)
}

CP_steps_fit <- function()
{
  OPTinf <- cleversearch(BICinf, lower=c(-4, 1), upper=c(0, 5),
                         ngrid=5, logscale=TRUE, verbose=FALSE)
  ## estimating mortality with optimal lambdas
  FITinf <- PSinfant(Y=Y1, E=E1, lambdas=OPTinf$par, WEI=WEI1)
  ## extract estimated linear predictor, log-mortality
  ETA1.hatI <- FITinf$ETA
  deltas <- deltasFUN(FITinf)
  return(list(deltas = deltas, OPTinf = OPTinf))
}
CP_steps_for <- function(obj)
{  
  deltas <- obj[["deltas"]]
  OPTinf <- obj[["OPTinf"]]
  ## where to apply the constraints
  S <- matrix(1, m, n)
  S[,1:n1] <- 0
  
  ## modelling with CP-splines
  FITcon <- CPSfunction(Y=Y_CP, E=E_CP, lambdas=OPTinf$par,
                        WEI=WEI, deltas=deltas, S=S,
                        verbose=FALSE)
  ## estimated and forecast linear predictor, log-mortality
  ETA.hatC <- FITcon$ETA
  
  rates_pred <- tail(x = ETA.hatC, n = c(m,h_step))
}

# Mortality data
load(here('output','mortality.Rdata'))

country <- "ita"
gender <- "woman"

Y <- eval(parse(text = paste('Y', country, gender, sep = '_')))
N <- eval(parse(text = paste('N', country, gender, sep = '_')))
Tmax <- nrow(Y)
Z <- ncol(Y)

n_for <- 25  
train <- 58  
h_step <- h <- 10 

method <- "Nelder-Mead"
parallel <- FALSE
maxcl <- 1
rep <- 1

years_bak <- years

LC <- lc(link = 'log')
APC <- apc(link = 'log')  
CBD <- cbd(link = 'log')
f2 <- function(x, ages) mean(ages) - x
constPlat <- function(ax, bx, kt, b0x, gc, wxt, ages){
  nYears <- dim(wxt)[2]
  x <- ages
  t <- 1:nYears
  c <- (1 - tail(ages, 1)):(nYears - ages[1])
  xbar <- mean(x)
  phiReg <- lm(gc ~ 1 + c + I(c ^ 2), na.action = na.omit)
  phi <- coef(phiReg)
  gc<-gc-phi[1]-phi[2]*c-phi[3]*c^2
  kt[2,]<-kt[2,]+2*phi[3]*t
  kt[1,]<-kt[1,]+phi[2]*t+phi[3]*(t^2-2*xbar*t)
  ax<-ax+phi[1]-phi[2]*x+phi[3]*x^2
  ci <- rowMeans(kt, na.rm = TRUE)
  ax<-ax+ci[1]+ci[2]*(xbar-x)
  kt[1, ] <- kt[1, ] - ci[1]
  kt[2, ] <- kt[2, ] - ci[2]
  return(list(ax=ax,bx=bx,kt=kt,b0x=b0x,gc=gc))
}
PLAT <- StMoMo(link = "log", staticAgeFun = TRUE,
               periodAgeFun = c("1", f2), cohortAgeFun = "1", constFun = constPlat)
RH <- rh(link = "log", cohortAgeFun = "1")


# Creating increasing sample size datasets
data_list <- lapply(train:Tmax,
                    function(t){
                      Y[1:t,]/N[1:t,]
                    })

Y_list <- lapply(train:Tmax,
                 function(t){
                   Y[1:t,]
                 })
N_list <- lapply(train:Tmax,
                 function(t){
                   N[1:t,]
                 })

times <- vector("list", length(data_list))
for(it in seq_along(data_list)){
  years <- years_bak
  ### Setup
  data <- data_list[[it]]
  Y <- Y_list[[it]]
  Y1 <- t(Y)
  N <- N_list[[it]]
  E1 <- t(Y)
  years <- years_it <- years1 <- years[1:nrow(Y)]
  ## BSP
  model <- bsp.model(data,
                     delta = delta,
                     age_knots = age_knots,
                     kernel = matern_kernel)
  ## LC
  stmomo_data <- structure(list(Dxt = t(Y), 
                                Ext = t(N), 
                                ages = ages, 
                                years = years_it, 
                                type = 'central', 
                                series = 'female',  
                                label = country), 
                           class = "StMoMoData")
  wxt_osa <- genWeightMat(ages = stmomo_data$ages,
                          years = years_it, clip = 3)
  ## CP
  years_CP <- c(years1, tail(years1,1) + 1:h)
  n <- length(years_CP)
  m <- length(ages)
  n1 <- length(years1)
  ETA1 <- log(Y1/E1)
  Y_CP <- matrix(10, m, n)
  Y_CP[1:m,1:n1] <- Y1
  E_CP <- matrix(10, m, n)
  E_CP[1:m,1:n1] <- E1
  WEI1 <- matrix(1, m, n1)
  WEI1[E1==0] <- 0
  WEI <- cbind(WEI1, matrix(0, m, n-n1))
  BICinf <- function(par){
    FITinf <- PSinfant(Y=Y1, E=E1, lambdas=par, WEI=WEI1)
    return(FITinf$bic)
  }
  
  mark(
    # BSP fit
    fit <- bsp.fit(model, 
                   rep = rep, 
                   method = method, 
                   parallel = parallel),
    # BSP forecast
    forecast <- forecast_v2(),
    # LC fit
    LC_fit <- fit(object = LC,
                  data = stmomo_data,
                  wxt = wxt_osa,
                  gc.order = NULL,
                  years = years_it,
                  verbose = FALSE),
    # LC forecast
    LC_for <- simulate(LC_fit, h = h, nsim = 100),
    # APC fit
    APC_fit <- fit(object = APC,
                   data = stmomo_data,
                   wxt = wxt_osa,
                   gc.order = c(1, 1, 0),
                   years = years_it,
                   verbose = FALSE),
    # APC forecast
    APC_for <- simulate(APC_fit, h = h, nsim = 100),
    # CBD fit
    CBD_fit <- fit(object = CBD,
                   data = stmomo_data,
                   wxt = wxt_osa,
                   gc.order = NULL,
                   years = years_it,
                   verbose = FALSE),
    # CBD forecast
    CBD_for <- simulate(CBD_fit, h = h, nsim = 100),
    # PLAT fit
    PLAT_fit <- fit(object = PLAT,
                   data = stmomo_data,
                   wxt = wxt_osa,
                   gc.order = c(2, 0, 0),
                   years = years_it,
                   verbose = FALSE),
    # PLAT forecast
    PLAT_for <- simulate(PLAT_fit, h = h, nsim = 100),
    # RH fit
    RH_fit <- fit(object = RH,
                  data = stmomo_data,
                  wxt = wxt_osa,
                  gc.order = c(2, 0, 0),
                  years = years_it,
                  verbose = FALSE),
    # RH forecast
    RH_for <- simulate(RH_fit, h = h, nsim = 100),
    # CP fit
    CP_fit <- CP_steps_fit(),
    # CP forecast
    CP_for <- CP_steps_for(CP_fit),
    # # HU fit
    HU_fit <- model.HyndmanUllah(data = t(data),
                                 x = ages,
                                 y = years_it),
    # # HU forecast
    HU_for <- predict(HU_fit, h = h_step, level = 95),
    check = FALSE,
    iterations = 10) -> tmp
  
  print("Done")
  
  times_it <- mutate(tmp, 
                     n = nrow(data),
                     label = c("BSP_fit", "BSP_for",
                               "LC_fit", "LC_for",
                               "APC_fit", "APC_for",
                               "CBD_fit", "CBD_for",
                               "PLAT_fit", "PLAT_for",
                               "RH_fit", "RH_for",
                               "CP_fit", "CP_for",
                               "HU_fit", "HU_for"))
  
  times[[it]] <- times_it
}

output_collector <- append(output_collector, 
                           list(warnings = warnings()))

# save(list = c('output_collector',
#               'times',
#               'data_list'),
#      file = here('output','comp_times.Rdata'))

times %>%
  map_dfr(.f = . %>%
            select(-expression) %>%
            relocate(label, .before = min) %>%
            select(label, median, `itr/sec`, n)) -> res

res %>%
  mutate(model = str_extract(label, pattern = "[A-Z]+"),
         what = str_extract(label, pattern = "[a-z]+$")) %>%
  ggplot(aes(x = n, y = median, color = model)) +
  geom_line() +
  facet_grid(row = vars(what))


