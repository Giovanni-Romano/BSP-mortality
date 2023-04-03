# Kernel function defining the correlation between the spline weights
matern_kernel <- function(x, phi = 0.5, kappa = 2){
  uphi <- x/phi
  uphi <- ifelse(x > 0,
                 (((2^(-(kappa-1)))/ifelse(0, Inf, base::gamma(kappa))) *
                    (uphi^kappa) *
                    besselK(x=uphi, nu=kappa)), 1)
  uphi[x > 600*phi] <- 0
  return(uphi)
}

# Delta value indicating frequency of observations (here equally spaced)
delta <- 0.01

# Spline knots
age_knots <- c(seq(from=5, to=40, by=5), 50, 60, seq(from=70, to=95, by=5))

# How many times to repeat MLE estimation
rep <- 5
