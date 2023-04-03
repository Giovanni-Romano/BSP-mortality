require(here)
require(tidyverse)
require(StMoMo)
require(doParallel)

rm(list=ls())

set.seed(5424)

load(here('output','mortality_extra.Rdata'))
output_collector <- list()
output_collector <- append(output_collector, list(info=Sys.time()))

# PLAT: reduced Plat model via StMoMo package
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

# Set up minimum training window, how many step-ahead to predict
train <- 50
h_step <- 10
output_collector <- append(output_collector, list(train = train, h_step = h_step))

# Function to fit  and forecast
fit_sim <- function(t, train, mod, h, stmomo_data, years, gc.order = NULL, nsim = 10){
  wxt_osa <- genWeightMat(ages = stmomo_data$ages,
                          years = years[1:t], clip = 3)
  
  err2 <- err <- try(log('a'), silent = TRUE) # random error
  while(inherits(err, "try-error") | inherits(err2, "try-error")){ # keep repeating if there is an error
    err <- try(fit <- fit(mod,
                          data = stmomo_data,
                          wxt = wxt_osa,
                          gc.order = gc.order,
                          years = years,
                          years.fit = years[1:t],
                          verbose = FALSE))
    print('fit done')
    err2 <- try(sim <- simulate(fit, h = h, nsim = nsim))
  }
  return(sim)
}

# Function to set up the rolling window: minimum window given by train, then 
# increasing the training set by 1 each time
rolling <- function(cg){
  country <- sub("_.*", "", cg)
  gender <- sub(".*_", "", cg)
  print(paste('Doing', country, gender))
  Y <- eval(parse(text = paste('Y', country, gender, sep = '_')))
  N <- eval(parse(text = paste('N', country, gender, sep = '_')))
  
  years <- years[1:nrow(Y)]
  
  stmomo_data <- structure(list(Dxt = t(Y), 
                                Ext = t(N), 
                                ages = ages, 
                                years = years, 
                                type = 'central', 
                                series = ifelse(gender == 'man', 'male', 'female'),  # whether "male", "female" or "total"
                                label = country), 
                           class = "StMoMoData")
  
  wxt <- genWeightMat(ages = stmomo_data$ages, 
                      years = stmomo_data$years, 
                      clip = 3)
  
  registerDoParallel(cores = 50)
  cl <- makeCluster(50, type = "FORK")
  PLATsim_ <- parLapply(cl,
                        train:length(years),
                        . %>%
                          fit_sim(., 
                                  train = train, 
                                  mod = PLAT, 
                                  h = h_step,
                                  stmomo_data = stmomo_data,
                                  years = years,
                                  gc.order = c(2, 0, 0),
                                  nsim = 100))
  stopCluster(cl)
  
  
  PLATsim <- PLATsim_ %>%
    modify(. %>% pluck('rates'))
  return(PLATsim)
}

# Set up parameter with country-gender combination
cg <- expand.grid(country = c('fra','dnk','cze'),
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
     file = here('output','PLAT_for_extra.Rdata'))


