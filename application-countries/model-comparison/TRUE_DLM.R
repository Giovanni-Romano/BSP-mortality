pred_pseudo <- predict(rwd_pseudo,
                       interval = 'prediction',
                       newdata = SSModel(matrix(NA, nrow = h, ncol = Z) ~ -1 +
                                           SSMcustom(Z = rwd_pseudo$Z, T = rwd_pseudo$T,
                                                     R = rwd_pseudo$R, Q = rwd_pseudo$Q,
                                                     state_names = rownames(rwd_pseudo$a1),
                                                     n = h), # total number of years, dimensionality check
                                         distribution = 'gaussian',
                                         H = rwd_pseudo$H),
                       level = 0.95,
                       type = 'response')
pseudo_smooth <- KFS(rwd_pseudo, smoothing = 'state')
V_rwd <- pseudo_smooth$V[1:(K+1),1:(K+1),1]
V_rwd_drift <- pseudo_smooth$V[-c(1:(K+1)),-c(1:(K+1)),1]
Qt_pred <- rwd_pseudo$Q[1:(K+1),1:(K+1),1]#V_rwd
H_rwd <- rwd_pseudo$H[,,1]
DELTA_var <- rwd_pseudo$Q[-c(1:(K+1)),-c(1:(K+1)),1]#dU_var_mean*Tt_pred^2 #diag(Tt_pred^2 * drift_median_var)
Tmat <- rwd_pseudo$T
a1_pred <- tail(pseudo_smooth$alphahat, n = 1 ,k = 1) # Now a1 is last observed year
U_pred <- matrix(NA, nrow = h, ncol = K+1)
varU_pred <- array(NA, dim = c(2*(K+1), 2*(K+1), h))
f_pred <- matrix(NA, nrow = h, ncol = Z)
varf_pred <- array(NA, dim = c(Z, Z, h))
for(tt in 1:h){
  if(tt == 1){
    U_pred[1,] <- as.numeric(a1_pred[1:(K+1)]) + Tt_pred%*%dU
    varU_pred[,,1] <- Tmat[,,1]%*%pseudo_smooth$V[,,1]%*%t(Tmat[,,1]) + rwd_pseudo$Q[,,1] #+ U_var_mean
    f_pred[1,] <- Zt_pred_%*%U_pred[1,]
    varf_pred[,,1] <- H_rwd + 
      rwd_pseudo$Z[,,1] %*% varU_pred[,,1] %*% t(rwd_pseudo$Z[,,1])
  } else {
    U_pred[tt,] <- U_pred[tt-1,] + Tt_pred%*%dU
    varU_pred[,,tt] <- Tmat[,,1]%*%varU_pred[,,tt-1]%*%t(Tmat[,,1]) + rwd_pseudo$Q[,,1]#+ U_var_mean
    f_pred[tt,] <- Zt_pred_%*%U_pred[tt,]
    varf_pred[,,tt] <- H_rwd + 
      rwd_pseudo$Z[,,1] %*% varU_pred[,,tt] %*% t(rwd_pseudo$Z[,,1])
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
