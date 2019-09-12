library(Rcpp)
library(RcppArmadillo)
sourceCpp("LinkedEmu_Cpp.cpp")

######Functions for prediction from a linked emulator using simulation/Monte Carlo method

#Sample handling
multimodel_samples <- function(nvar, predpts, loc_vec, y1_samp, designs, b_vecs, res_vecs, F_mats, corr_mats, cov_type, scale_params){
  
  sampsizes <- length(y1_samp)
  samps <- matrix(nrow=sampsizes, ncol=nvar)
  samps[,1] <- y1_samp
  for(i in 1:(nvar-1)){
    if(!is.null(scale_params)){
      samps[,i] <- (samps[,i] - scale_params[[i]][1]) / scale_params[[i]][2]
    }
    if(predpts[i]=="N"){
      samps[,i+1] <- GPsamp_multiple(1, as.matrix(samps[,i]), designs[[i+1]], b_vecs[[i+1]], res_vecs[[i+1]], F_mats[[i+1]], corr_mats[[i+1]], cov_type)
    }
    else{
      samps[,i+1] <- condsamps(samps[,i], loc_vec[i], designs[[i+1]], b_vecs[[i+1]], res_vecs[[i+1]], F_mats[[i+1]], corr_mats[[i+1]], as.numeric(predpts[i]), cov_type)
    }
  }
  return(samps)
}

#Prediction
MM_pred_simu <- function(sampsize, predpts, loc_vec, designs, b_vecs, res_vecs, F_mats, corr_mats, cov_type=c("Gaussian", "Matern_32", "Matern_52"), scale_params=NULL){
  
  for(i in 1:length(designs)){
    designs[[i]] <- as.matrix(designs[[i]])
    if(i <= length(predpts)){
      if(!is.null(predpts[[i]])){
        predpts[[i]] <- as.matrix(predpts[[i]]) 
      }
    }
  }
  
  npts <- nrow(predpts[[1]])
  num_models <- length(loc_vec)+1
  GP_sample <- matrix(nrow=npts, ncol=sampsize)
  
  pred_x1 <- predpts[[1]]
  xn_1 <- designs[[1]]
  b1 <- b_vecs[[1]]
  yn_1 <- res_vecs[[1]]
  F_mat1 <- F_mats[[1]]
  corr_mat1 <- corr_mats[[1]]
  
  cov_type <- match.arg(cov_type)
  cov_int <- 1
  if(cov_type=="Matern_32"){
    cov_int <- 2
  }
  else if(cov_type=="Matern_52"){
    cov_int <- 3
  }
  
  y1_samps_all <- GPsamp_multiple(sampsize, pred_x1, xn_1, b1, yn_1, F_mat1, corr_mat1, cov_int)
  
  for(i in 1:npts){
    if(length(predpts) > 1){
      if(!is.null(predpts[[2]])){
        predpts_curr <- t(as.matrix(predpts[[2]])[i,])
      }
      else{
        predpts_curr <- "N"
      }
    }
    else{
      predpts_curr <- "N"
    }
    if(num_models>2){
      for(j in 2:(num_models-1)){
        if(j <= length(predpts)){
          if(!is.null(predpts[[j+1]])){
            predpts_curr <- rbind(predpts_curr, as.matrix(predpts[[j+1]])[i,])
          }
          else{
            predpts_curr <- rbind(predpts_curr, "N")
          }
        }
        else{
          predpts_curr <- rbind(predpts_curr, "N")
        }
      }
    }
    y1_samp <- y1_samps_all[i,]
    samps_all <- multimodel_samples(num_models, predpts_curr, loc_vec, y1_samp, designs, b_vecs, res_vecs, F_mats, corr_mats, cov_int, scale_params)
    samp_final <- samps_all[, num_models]
    GP_sample[i,] <- samp_final
  }
  
  return(GP_sample)
}

######Functions for prediction from a linked emulator using theoretical method

#Mean and variance of uncertain input
MM_theor_setup <- function(xnew, xn_1, y1, beta0_1, inv_corrmat1, b1, sig2_1, nugget_1, scale_params){
  
  xn_1 <- as.matrix(xn_1)
  npts <- nrow(inv_corrmat1)
  nin1 <- length(b1)
  if(nin1<length(xnew)){
    xnew_1 <- xnew[1:nin1]
    xnew_2 <- (xnew[-(1:nin1)])
    xnew_2 <- as.vector(xnew_2)
  }
  else{
    xnew_1 <- xnew
    xnew_2 <- NULL
  }
  c1 <- vector(length=npts)
  for(i in 1:npts){
    c1[i] <- corrfunc_Gauss_Cpp(b1, xn_1[i,] - xnew_1)
  }
  mu1 <- beta0_1 + c1 %*% (inv_corrmat1 %*% (y1 - beta0_1))
  sig2_y1 <- sig2_1 * (1 + nugget_1 - c1 %*% inv_corrmat1 %*% c1)
  
  if(!is.null(scale_params)){
    mu1 <- (mu1 - scale_params[1]) / scale_params[2]
    sig2_y1 <- sig2_y1 / scale_params[2]^2
  }
  
  return(list(xnew_2=xnew_2, mu1=mu1, sig2_y1=sig2_y1))
}

#Mean and variance of linked emulator for two-model chain
PVar_full <- function(predpts, xn_1, xn_2, y1, y2, beta0_1, beta0_2, inv_corrmat1, inv_corrmat2, b1, b2, sig2_1, sig2_2, nugget_1, nugget_2, scale_params=NULL){
  
  xn_1 <- as.matrix(xn_1)
  xn_2 <- as.matrix(xn_2)
  predpts <- as.matrix(predpts)
  a <- inv_corrmat2 %*% (y2 - beta0_2)
  npts <- nrow(predpts)
  Yn_Exps <- vector(length=npts)
  Yn_Vars <- vector(length=npts)
  
  for(i in 1:npts){
    xnew <- predpts[i,]
    setup_vals <- MM_theor_setup(xnew, xn_1, y1, beta0_1, inv_corrmat1, b1, sig2_1, nugget_1, scale_params)
    xnew_2 <- setup_vals$xnew_2
    mu1 <- setup_vals$mu1
    sig2_y1 <- setup_vals$sig2_y1
    
    if(is.null(xnew_2)){
      Yn_Exps[i] <- PExp_mult_Cpp_nonewin(xn_2, a, beta0_2, b2, mu1, sig2_y1) 
      Yn_Vars[i] <- PVar_mult_Cpp_nonewin(xn_2, a, b2, inv_corrmat2, mu1, sig2_y1, sig2_2, nugget_2)
    }
    else{
      Yn_Exps[i] <- PExp_mult_Cpp_full(xn_2, a, beta0_2, b2, mu1, sig2_y1, xnew_2) 
      Yn_Vars[i] <- PVar_mult_Cpp_full(xn_2, a, b2, inv_corrmat2, mu1, sig2_y1, sig2_2, nugget_2, xnew_2)
    }
  }
  
  return(list(means=Yn_Exps, vars=Yn_Vars))
}

#Prediction from a chain of any length
MM_pred_theory <- function(predpts, designs, b_vecs, res_vecs, inv_corrmats, beta0s, sig2s, nuggets, scale_params=NULL){
  
  num_models <- length(designs)
  
  for(i in 1:num_models){
    designs[[i]] <- as.matrix(designs[[i]])
    if(i <= length(predpts)){
      if(!is.null(predpts[[i]])){
        predpts[[i]] <- as.matrix(predpts[[i]]) 
      }
    }
  }
  
  npts <- nrow(predpts[[1]])
  
  pred_x1 <- predpts[[1]]
  xn_1 <- designs[[1]]
  b1 <- b_vecs[[1]]
  y1 <- res_vecs[[1]]
  inv_corrmat1 <- inv_corrmats[[1]]
  beta0_1 <- beta0s[[1]]
  sig2_1 <- sig2s[[1]]
  nugget_1 <- nuggets[[1]]
  
  xn_2 <- designs[[2]]
  b2 <- b_vecs[[2]]
  y2 <- res_vecs[[2]]
  inv_corrmat2 <- inv_corrmats[[2]]
  beta0_2 <- beta0s[[2]]
  sig2_2 <- sig2s[[2]]
  nugget_2 <- nuggets[[2]]
  
  if(length(predpts) > 1){
    if(!is.null(predpts[[2]])){
      pred_x2 <- predpts[[2]]
      predpts_m12 <- cbind(pred_x1, pred_x2)
    }
    else{
      predpts_m12 <- pred_x1
    }
  }
  else{
    predpts_m12 <- pred_x1
  }
  
  if(is.null(scale_params)){
    res1 <- PVar_full(predpts_m12, xn_1, xn_2, y1, y2, beta0_1, beta0_2, inv_corrmat1, inv_corrmat2, b1, b2, sig2_1, sig2_2, nugget_1, nugget_2)
  }
  else{
    res1 <- PVar_full(predpts_m12, xn_1, xn_2, y1, y2, beta0_1, beta0_2, inv_corrmat1, inv_corrmat2, b1, b2, sig2_1, sig2_2, nugget_1, nugget_2, scale_params[[1]])
  }
  
  E_curr <- res1$mean
  Var_curr <- res1$var
  
  if(num_models>2){
    for(j in 3:(num_models)){
      
      xn_new <- designs[[j]]
      bnew <- b_vecs[[j]]
      ynew <- res_vecs[[j]]
      inv_corrmatnew <- inv_corrmats[[j]]
      beta0_new <- beta0s[[j]]
      sig2_new <- sig2s[[j]]
      nugget_new <- nuggets[[j]]
      anew <- inv_corrmatnew %*% (ynew - beta0_new)
      
      for(i in 1:npts){
        if(j <= length(predpts)){
          if(!is.null(predpts[[j]])){
            pred_new <- predpts[[j]][i,]
          }
          else{
            pred_new <- "N"
          }
        }
        else{
          pred_new <- "N"
        }
        if(pred_new=="N"){
          mu_y <- E_curr[i]
          sig2_y <- Var_curr[i]
          if(!is.null(scale_params)){
            mu_y <- (mu_y - scale_params[[j-1]][1]) / scale_params[[j-1]][2]
            sig2_y <- sig2_y / scale_params[[j-1]][2]^2
          }
          E_curr[i] <- PExp_mult_Cpp_nonewin(xn_new, anew, beta0_new, bnew, mu_y, sig2_y) 
          Var_curr[i] <- PVar_mult_Cpp_nonewin(xn_new, anew, bnew, inv_corrmatnew, mu_y, sig2_y, sig2_new, nugget_new)
        }
        else{
          
          E_curr[i] <- PExp_mult_Cpp_full(xn_new, anew, beta0_new, bnew, E_curr[i], Var_curr[i], pred_new) 
          Var_curr[i] <- PVar_mult_Cpp_full(xn_new, anew, bnew, inv_corrmatnew, E_curr[i], Var_curr[i], sig2_new, nugget_new, pred_new)
        }
      }
    }
  }
  Var_curr[which(Var_curr<0)] <- 0
  return(list(means=E_curr, vars=Var_curr))
}

######Functions for sensitivity analysis - one emulator

#Posterior mean of the expectation with no fixed inputs
PExp_EY <- function(b, xn, beta_hat, inv_corrmat, F_mat, yn, sample, is.exact=T, G=NULL, S=NULL){
  
  if(is.exact==T){
    R_mat <- Rint_dG_Cpp_exact(as.matrix(sample))
    T_mat <- Tint_dG_Cpp_exact(as.matrix(sample), b, xn)
  }
  else{
    R_mat <- Rint_dG_Cpp_importance(G, S, as.matrix(sample))
    T_mat <- Tint_dG_Cpp_importance(G, S, as.matrix(sample), b, xn)
  }
  
  e <- inv_corrmat %*% (yn - F_mat * beta_hat)
  return(R_mat %*% beta_hat + T_mat %*% e)
}

#Posterior mean of the expectation given some fixed inputs
PExp_Egxp <- function(b, xn, xp, xp_loc, beta_hat, e, sample, is.exact=T, Gdp=NULL, Sdp=NULL){
  
  if(is.exact==T){
    Rp_mat <- Rint_xp_Cpp_exact(as.matrix(sample), xp, xp_loc)
    Tp_mat <- Tint_xp_Cpp_exact(as.matrix(sample), b, xn, xp, xp_loc)
  }
  else{
    Rp_mat <- Rint_xp_Cpp_importance(Gdp, Sdp, as.matrix(sample), xp, xp_loc)
    Tp_mat <- Tint_xp_Cpp_importance(Gdp, Sdp, as.matrix(sample), b, xn, xp, xp_loc)
  }
  return(Rp_mat %*% beta_hat + Tp_mat %*% e)
}

#Posterior variance of the expectation given some fixed inputs
PVar_Egxi <- function(inv_Cmat, F_mat, b, xn, W, nugget, xi, xi_loc, sig2hat, sample1, sample2, is.exact=T, Gdi=NULL, Sdi=NULL){

  if(is.exact==T){
    Cstar_term <- Cstarint_xp2d_Cpp_exact(as.matrix(sample1), as.matrix(sample2), inv_Cmat, F_mat, b, xn, W, nugget, xi, xi_loc)
  }
  else{
    Cstar_term <- Cstarint_xp2d_Cpp_importance(Gdi, Sdi, as.matrix(sample1), as.matrix(sample2), inv_Cmat, F_mat, b, xn, W, nugget, xi, xi_loc)
  }
  return(sig2hat * Cstar_term)
}

#Posterior variance of the expectation with no fixed inputs
PVar_EY <- function(inv_Cmat, F_mat, b, xn, W, nugget, sig2hat, sample1, sample2, is.exact=T, G=NULL, S=NULL){
  
  if(is.exact==T){
    Cstar_term <- Cstarint_dGdG_Cpp_exact(as.matrix(sample1), as.matrix(sample2), inv_Cmat, F_mat, b, xn, W, nugget)
  }
  else{
    Cstar_term <- Cstarint_dGdG_Cpp_importance(G, S, as.matrix(sample1), as.matrix(sample2), inv_Cmat, F_mat, b, xn, W, nugget)
  }
  return(sig2hat * Cstar_term)
}

#Posterior mean of the square of the expectation with no fixed inputs
PExp_EYsq <- function(inv_Cmat, F_mat, b, xn, W, nugget, e, beta_hat, sig2hat, sample, is.exact=T, G=NULL, S=NULL){
  if(is.exact==T){
    term1a <- Cstarint_dG_Cpp_exact(as.matrix(sample), inv_Cmat, F_mat, b, xn, W, nugget)
    term2 <- Mstarint_dG_Cpp_exact(as.matrix(sample), b, xn, e, beta_hat)
  }
  else{
    term1a <- Cstarint_dG_Cpp_importance(G, S, as.matrix(sample), inv_Cmat, F_mat, b, xn, W, nugget)
    term2 <- Mstarint_dG_Cpp_importance(as.matrix(sample), b, xn, e, beta_hat, G, S)
  }
  term1 <- sig2hat * term1a
  return(term1 + term2)
}

#Posterior mean of total variance
PExp_VarY <- function(inv_Cmat, F_mat, b, xn, W, nugget, e, beta_hat, sig2hat, EsEY, VarsEY, sample, is.exact=T, G=NULL, S=NULL){
  if(is.exact==T){
    term1 <- PExp_EYsq(inv_Cmat, F_mat, b, xn, W, nugget, e, beta_hat, sig2hat, sample)
  }
  else{
    term1 <- PExp_EYsq(inv_Cmat, F_mat, b, xn, W, nugget, e, beta_hat, sig2hat, sample, is.exact=F, G=G, S=S)
  }
  res <- term1 - EsEY^2 - VarsEY
  return(res)
}

#Posterior mean of the expectation of square of the expectation with some inputs fixed
PExp_EEygxpsq <- function(inv_Cmat, F_mat, b, xn, W, nugget, e, beta_hat, sig2hat, xp_loc, sample1, sample2, sample3, is.exact=T, Gdp=NULL, Gp=NULL, Sdp=NULL, Sp=NULL){
  if(is.exact==T){
    term1a <- Cstarint_xp3d_Cpp_exact(as.matrix(sample1), as.matrix(sample2), as.matrix(sample3), inv_Cmat, F_mat, b, xn, W, nugget, xp_loc)
    term2 <- Mstarint_xp_Cpp_exact(as.matrix(sample1), as.matrix(sample2), as.matrix(sample3), b, xn, e, beta_hat, xp_loc)
  }
  else{
    term1a <- Cstarint_xp3d_Cpp_importance(Gdp, Sdp, Gp, Sp, as.matrix(sample1), as.matrix(sample2), as.matrix(sample3), inv_Cmat, F_mat, b, xn, W, nugget, xp_loc)
    term2 <- Mstarint_xp_Cpp_importance(Gdp, Sdp, Gp, Sp, as.matrix(sample1), as.matrix(sample2), as.matrix(sample3), b, xn, e, beta_hat, xp_loc)
  }
  term1 <- sig2hat * term1a
  return(term1 + term2)
}

#Posterior mean of Sobol' index
PExp_Vp <- function(inv_Cmat, F_mat, b, xn, W, nugget, e, beta_hat, sig2hat, EsEY, VarsEY, xp_loc, sample1, sample2, sample3, is.exact=T, Gdp=NULL, Gp=NULL, Sdp=NULL, Sp=NULL){
  if(is.exact==T){
    term1 <- PExp_EEygxpsq(inv_Cmat, F_mat, b, xn, W, nugget, e, beta_hat, sig2hat, xp_loc, as.matrix(sample1), as.matrix(sample2), as.matrix(sample3))
  }
  else{
    term1 <- PExp_EEygxpsq(inv_Cmat, F_mat, b, xn, W, nugget, e, beta_hat, sig2hat, xp_loc, as.matrix(sample1), as.matrix(sample2), as.matrix(sample3), is.exact=F, Gdp=Gdp, Gp=Gp, Sdp=Sdp, Sp=Sp)
  }
  res <- term1 - EsEY^2 - VarsEY
  return(res)
}

#Sobol' indices for all inputs
CalcSi <- function(inv_corrmat, F_mat, b, xn, W, nugget, e, beta_hat, sig2_hat, EsEY, VarsEY, EsVarY, sample1, sample2, samp_complete=T, sample3=NULL){
  nvar <- length(b)
  res <- vector(length=nvar)
  for(i in 1:(nvar)){
    if(samp_complete==T){
      samp1_curr <- sample1[,-i]
      samp2_curr <- sample2[,-i]
      samp3_curr <- as.matrix(sample1[,i])
    }
    else{
      samp1_curr <- sample1
      samp2_curr <- sample2
      samp3_curr <- sample3
    }
    res[i] <- PExp_Vp(inv_corrmat, F_mat, b, xn, W, nugget, e, beta_hat, sig2_hat, EsEY, VarsEY, i, samp1_curr, samp2_curr, samp3_curr)
    if (res[i] < 0){
      res[i] <- 0
    }
  }
  Si_hat_all <- res / rep(EsVarY, nvar)
  return(Si_hat_all)
}

#Sobol' indices for all two-way interactions
CalcSij <- function(inv_corrmat, F_mat, b, xn, W, nugget, e, beta_hat, sig2_hat, EsEY, VarsEY, EsVarY, Si, sample1, sample2, samp_complete=T, sample3=NULL){
  nvar <- length(b)
  Sij_hat_all <- matrix(nrow=nvar-1, ncol=nvar)
  for(i in 1:(nvar-1)){
    for(j in 1:(nvar)){
      if(i < j){
        if(samp_complete==T){
          samp1_curr <- sample1[,c(-i, -j)]
          samp2_curr <- sample2[,c(-i, -j)]
          samp3_curr <- as.matrix(sample1[,c(i, j)])
        }
        else{
          samp1_curr <- sample1
          samp2_curr <- sample2
          samp3_curr <- sample3
        }
        term1a <- PExp_Vp(inv_corrmat, F_mat, b, xn, W, nugget, e, beta_hat, sig2_hat, EsEY, VarsEY, c(i, j), samp1_curr, samp2_curr, samp3_curr)
        term1 <- term1a / EsVarY
        Sij_hat_all[i, j] <- term1 - Si[i] - Si[j]
        if(Sij_hat_all[i, j] < 0){
          Sij_hat_all[i, j] <- 0
        }
      }
    }
  }
  return(Sij_hat_all)
}

#Calculate and plot posterior means across the range of all inputs
CalcME <- function(b, xn, beta_hat, e, sample, seq_length, samp_complete=T, lim_low=NULL, lim_up=NULL, varnames=NULL){
  
  nvar <- length(b)
  name_flag <- 0
  if(is.null(varnames)){
    varnames <- vector(length=nvar)
    name_flag <- 1
  }
  PM_EYxi <- matrix(nrow=nvar, ncol=seq_length)
  xi_seq <- matrix(nrow=nvar, ncol=seq_length)
  if(is.null(lim_low)){
    lim_low <- rep(0, nvar)
  }
  if(is.null(lim_up)){
    lim_up <- rep(1, nvar)
  }
  for(i in 1:(nvar)){
    xi_seq[i,] <- seq(lim_low[i], lim_up[i], len=seq_length)
    if(samp_complete==T){
      samp_curr <- sample[,-i]
    }
    else{
      samp_curr <- sample
    }
    for(j in 1:seq_length){
      xi_curr <- xi_seq[i,j]
      PM_EYxi[i, j] <- PExp_Egxp(b, xn, xi_curr, i, beta_hat, e, samp_curr)
    }
  }
  
  plot(xi_seq[1,], PM_EYxi[1, ], type="l", xlab="xi", ylab="E(Y|xi)", ylim=c(min(PM_EYxi), max(PM_EYxi)), xlim=c(min(lim_low), max(lim_up)))
  if(name_flag==1){
    varnames[1] <- paste("x", 1, sep="")
  }
  
  if(nvar>1){
    for(i in 2:nvar){
      lines(xi_seq[i,], PM_EYxi[i, ], col=i)
      if(name_flag==1){
        varnames[i] <- paste("x", i, sep="")
      }
    }
  }
  legend("topright", legend=varnames, lty=1, col=1:nvar, cex=0.8)
  
  return(PM_EYxi)
}

#Plot posterior variance bounds of the mean given an input
PlotMEBounds <- function(varnum, varmean, inv_corrmat, F_mat, b, xn, W, nugget, beta_hat, e, sig2_hat, sample1, sample2, seq_length, lim_low=0, lim_up=1, samp_complete=T){
  xi_seq <- seq(lim_low, lim_up, len=seq_length)
  PVar_Eyxi <- vector(length=seq_length)
  
  if(samp_complete==T){
    samp_curr1 <- sample1[,-varnum]
    samp_curr2 <- sample2[,-varnum]
  }
  else{
    samp_curr1 <- sample1
    samp_curr2 <- sample2
  }
  
  for(j in 1:seq_length){
    xi_curr <- xi_seq[j]
    PVar_Eyxi[j] <- PVar_Egxi(inv_corrmat, F_mat, b, xn, W, nugget, xi_curr, varnum, sig2_hat, samp_curr1, samp_curr2)
  }
  
  Stdevs <- sqrt(PVar_Eyxi)
  bounds_upper <- varmean + 2 * Stdevs
  bounds_lower <- varmean - 2 * Stdevs
  
  plot(xi_seq, bounds_upper, type="l", xlab=paste("x", varnum, sep=""), ylab=paste("E(Y|x", varnum, ")", sep=""), ylim=c(min(bounds_lower), max(bounds_upper)))
  lines(xi_seq, bounds_lower)
}

#Plot posterior expectation of interaction between two inputs
PlotInt <- function(xp_loc, b, xn, beta_hat, e, sample1, seq_length, lim_low=0, lim_up=1, samp_complete=T){
  xi_seq <- seq(lim_low, lim_up, len=seq_length)
  PM_EYxixj <- matrix(nrow=seq_length, ncol=seq_length)
  
  if(samp_complete==T){
    samp_curr <- sample1[,-xp_loc]
  }
  else{
    samp_curr <- sample1
  }
  
  for(i in 1:seq_length){
    for(j in 1:seq_length){
      xp_curr <- c(xi_seq[i], xi_seq[j])
      PM_EYxixj[i, j] <- PExp_Egxp(b, xn, xp_curr, xp_loc, beta_hat, e, samp_curr)
    }
  }
  filled.contour(xi_seq, xi_seq, PM_EYxixj, xlab=paste("x", xp_loc[1], sep=""), ylab=paste("x", xp_loc[2], sep=""))
}

#Fully flexible sensitivity analysis for one emulator
SensGP <- function(design, output, inv_corrmat, F_mat, b, nugget, sample1, sample2, samp_complete=T, sample3=NULL, seq_length=200, bounds=NULL, ints=NULL, lim_low=NULL, lim_up=NULL, varnames=NULL){
  
  GP_quantities <- GP_Posterior_values_weak(design, output, inv_corrmat, F_mat)
  beta_hat <- GP_quantities$beta_hat
  sig2_hat <- GP_quantities$sig2_hat
  
  W_inv <- t(F_mat) %*% inv_corrmat %*% F_mat
  W <- chol2inv(chol(W_inv))
  R_mat <- Rint_dG_Cpp_exact(sample1)
  T_mat <- Tint_dG_Cpp_exact(sample1, b, design)
  e <- inv_corrmat %*% (output - F_mat * beta_hat)
  
  EsEY <- R_mat %*% beta_hat + T_mat %*% e
  VarsEY <- PVar_EY(inv_corrmat, F_mat, b, design, W, nugget, sig2_hat, sample1, sample2)
  EsVarY <- PExp_VarY(inv_corrmat, F_mat, b, design, W, nugget, e, beta_hat, sig2_hat, EsEY, VarsEY, sample1)
  
  ME <- CalcME(b, design, beta_hat, e, sample1, seq_length, samp_complete, lim_low, lim_up, varnames)
  
  Si_hat_all <- CalcSi(inv_corrmat, F_mat, b, design, W, nugget, e, beta_hat, sig2_hat, EsEY, VarsEY, EsVarY, sample1, sample2)
  Sij_hat_all <- CalcSij(inv_corrmat, F_mat, b, design, W, nugget, e, beta_hat, sig2_hat, EsEY, VarsEY, EsVarY, Si_hat_all, sample1, sample2)
  
  if(!is.null(bounds)){
    
    if(is.null(lim_low)){
      lim_low <- 0
    }
    if(is.null(lim_up)){
      lim_up <- 1
    }
    
    for(i in 1:length(bounds)){
      readline(prompt="Hit <return> to see next plot")
      var_curr <- bounds[i]
      E_curr <- ME[var_curr,]
      PlotMEBounds(var_curr, E_curr, inv_corrmat, F_mat, b, design, W, nugget, beta_hat, e, sig2_hat, sample1, sample2, seq_length, lim_low, lim_up, samp_complete)
    }
  }
  
  if(!is.null(ints)){
    
    if(is.null(lim_low)){
      lim_low <- 0
    }
    if(is.null(lim_up)){
      lim_up <- 1
    }
    
    if(length(ints)==2){
      readline(prompt="Hit <return> to see next plot")
      PlotInt(ints, b, design, beta_hat, e, sample1, seq_length, lim_low, lim_up, samp_complete)
    }
    else{
      for(i in 1:nrow(ints)){
        readline(prompt="Hit <return> to see next plot")
        PlotInt(ints[i,], b, design, beta_hat, e, sample1, seq_length, lim_low, lim_up, samp_complete)
      }
    }
  }  
  return(list(MainEffects=Si_hat_all, TwoWayInteractions=Sij_hat_all))
}

######Functions for sensitivity analysis for a linked emulator (theoretical method)

#Sensitivity analysis for the final model in a chain
SensFinal <- function(designs, b_vecs, res_vecs, inv_corrmats, beta0s, sig2s, nuggets, samples1, samples2, scale_params=NULL, seq_length=200, bounds=NULL, ints=NULL, lim_low=NULL, lim_up=NULL, varnames=NULL){
  
  nmod <- length(b_vecs)
  
  des_final <- designs[[nmod]]
  b_final <- b_vecs[[nmod]]
  res_final <- res_vecs[[nmod]]
  inv_corrmat_final <- inv_corrmats[[nmod]]
  beta0_final <- beta0s[[nmod]]
  sig2_final <- sig2s[[nmod]]
  nugget_final <- nuggets[[nmod]]
  samples1_final_xs <- samples1[[nmod]]
  samples2_final_xs <- samples2[[nmod]]
  F_mat_final <- as.matrix(rep(1, length(designs[[nmod]][,1])))
  
  if(nmod==2){
    F_mat1 <- as.matrix(rep(1, length(designs[[1]][,1])))
    corr_mat1 <- chol2inv(chol(inv_corrmats[[1]]))
    sample1_final_y <- GPsamp_multiple(1, samples1[[1]], designs[[1]], b_vecs[[1]], res_vecs[[1]], F_mat1, corr_mat1, 1)
    sample2_final_y <- GPsamp_multiple(1, samples2[[1]], designs[[1]], b_vecs[[1]], res_vecs[[1]], F_mat1, corr_mat1, 1)
    if(!is.null(scale_params)){
      sample1_final_y <- (sample1_final_y - scale_params[1]) / scale_params[2]
      sample2_final_y <- (sample2_final_y - scale_params[1]) / scale_params[2]
    }
  }
  else{
    designs[[nmod]] <- NULL
    b_vecs[[nmod]] <- NULL
    res_vecs[[nmod]] <- NULL
    inv_corrmats[[nmod]] <- NULL
    beta0s[[nmod]] <- NULL
    sig2s[[nmod]] <- NULL
    nuggets[[nmod]] <- NULL
    samples1[[nmod]] <- NULL
    samples2[[nmod]] <- NULL
    sample1_final_y <- MM_pred_theory(samples1, designs, b_vecs, res_vecs, inv_corrmats, beta0s, sig2s, nuggets, scale_params)
    sample2_final_y <- MM_pred_theory(samples2, designs, b_vecs, res_vecs, inv_corrmats, beta0s, sig2s, nuggets, scale_params)
  }
  
  sample1_final <- cbind(sample1_final_y, samples1_final_xs)
  sample2_final <- cbind(sample2_final_y, samples2_final_xs)
  
  res1 <- SensGP(des_final, res_final, inv_corrmat_final, F_mat_final, b_final, nugget_final, sample1_final, sample2_final, seq_length=seq_length, bounds=bounds, ints=ints, lim_low=lim_low, lim_up=lim_up, varnames=varnames)
  return(res1)
}

#Posterior mean of expectation with some controllable inputs to the chain fixed
Multimod_EsEygxp <- function(xp, xp_loc, sample, xn_1, xn_2, y1, y2, beta0_1, beta0_2, inv_corrmat1, inv_corrmat2, b1, b2, sig2_1, sig2_2, nugget_1, nugget_2, scale_params=NULL){
  
  nsamps <- nrow(sample)
  num_fixed <- length(xp_loc)
  means <- vector(length=nsamps)
  
  for(i in 1:nsamps){
    
    xnew <- sample[i,]
    for(j in 1:num_fixed){
      loc_curr <- xp_loc[j]
      if(loc_curr==1){
        xnew <- append(xp[j], xnew)
      }
      else{
        xnew <- append(xnew, xp[j], after=loc_curr-1)
      }
    }
    
    setup_vals <- MM_theor_setup(xnew, xn_1, y1, beta0_1, inv_corrmat1, b1, sig2_1, nugget_1, scale_params)
    xnew_2 <- setup_vals$xnew_2
    mu1 <- setup_vals$mu1
    sig2_y1 <- setup_vals$sig2_y1
    if(is.null(xnew_2)){
      means[i] <- PExp_mult_Cpp_nonewin(xn_2, y2, beta0_2, b2, inv_corrmat2, mu1, sig2_y1) 
    }
    else{
      means[i] <- PExp_mult_Cpp_full(xn_2, y2, beta0_2, b2, inv_corrmat2, mu1, sig2_y1, xnew_2)
    }
  }
  
  EsEygxp <- mean(means)
  return(EsEygxp)
}

#Calculate and plot posterior means across the range of all controllable inputs to a chain
Multimod_calcME <- function(sample, seq_length, xn_1, xn_2, y1, y2, beta0_1, beta0_2, inv_corrmat1, inv_corrmat2, b1, b2, sig2_1, sig2_2, nugget_1, nugget_2, samp_complete=T, lim_low=NULL, lim_up=NULL, varnames=NULL, scale_params=NULL){
  
  nvar <- ncol(sample)
  name_flag <- 0
  if(is.null(varnames)){
    varnames <- vector(length=nvar)
    name_flag <- 1
  }
  
  xi_seq <- matrix(nrow=nvar, ncol=seq_length)
  PM_Eyxi <- matrix(nrow=nvar, ncol=seq_length)
  PVar_Eyxi <- matrix(nrow=nvar, ncol=seq_length)
  if(is.null(lim_low)){
    lim_low <- rep(0, nvar)
  }
  if(is.null(lim_up)){
    lim_up <- rep(1, nvar)
  }
  for(i in 1:(nvar)){
    xi_seq[i,] <- seq(lim_low[i], lim_up[i], len=seq_length)
    if(samp_complete==T){
      samp_curr <- sample[,-i]
    }
    else{
      samp_curr <- sample
    }
    for(j in 1:seq_length){
      xi_curr <- xi_seq[i,j]
      res_temp <- Multimod_EsEygxp(xi_curr, i, samp_curr, xn_1, xn_2, y1, y2, beta0_1, beta0_2, inv_corrmat1, inv_corrmat2, b1, b2, sig2_1, sig2_2, nugget_1, nugget_2, scale_params)
      PM_Eyxi[i, j] <- res_temp
    }
  }
  
  plot(xi_seq[1,], PM_Eyxi[1, ], type="l", xlab="xi", ylab="E(Y|xi)", ylim=c(min(PM_Eyxi), max(PM_Eyxi)), xlim=c(min(lim_low), max(lim_up)))
  if(name_flag==1){
    varnames[1] <- paste("x", 1, sep="")
  }
  
  if(nvar>1){
    for(i in 2:nvar){
      lines(xi_seq[i,], PM_Eyxi[i, ], col=i)
      if(name_flag==1){
        varnames[i] <- paste("x", i, sep="")
      }
    }
  }
  legend("topright", legend=varnames, lty=1, col=1:nvar, cex=0.8)
  
  return(PM_Eyxi)
}
