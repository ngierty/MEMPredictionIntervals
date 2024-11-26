####################################################################################
# Functions for performing conformal prediction for the measurement error model
####################################################################################

library(MASS) # for generalized inverse
library(sn) # for vech function
library(quadprog) # for solving QP

####################################
# Functions for generating the data
####################################

calc_design_matrix <- function(x, p, zeta){
  
  n <- length(x)
  if(is.null(zeta)){
    K <- 0
  }else{
    K <- length(zeta)
  }
  X <- matrix(NA, nrow = n, ncol = p + K + 1)
  X[, 1] <- 1
  
  for(i in 2:(p+1)){
    if(p > 1){warning('functions to construct prediction intervals only made for p=1')}
    
    X[, i] <- x^(i-1)
  }
  
  if(K > 0){
    for(i in (p+2):(p+K+1)){
      j <- i - (p+1)
      X[, i] <- ((x - zeta[j])^p)*((x - zeta[j]) > 0)
    }
  }
  return(X)
}


#' Generate the data used to demonstrate that the proposed non-conformity
#' scores can be used in the measurement error model setting
#'
#' \code{gen_me_data} Generate measurement error model data
#'
#' @param nsamp number of observed and test samples
#' @param beta vector of regression parameters
#' @param sigma2 intrinsic scatter in response y; see equation (1) of manuscript
#' @param x_seed seed for true covariate values
#' @param data_seed seed for observed covariate and response values
#' @param delta_l lower bound for variance dependence on true covariates
#' @param delta_u upper bound for variance dependence on true covariates
#' @param dist distribution for the errors
#' @param known_vars flag for whether the variances are assumed to be known (up to multiple)
#' @return Returns the observed covariate, response, and standard errors; note "_train" combines "_train_not_cal" and "_cal"
#' @export
gen_me_data <- function(opt){
  
  # Generate the beta vector from the input
  beta <- c(as.double(strsplit(opt$beta, ',')[[1]]))
  
  # Generate the tau vector from the input
  tau_str <- strsplit(opt$tau, ',')[[1]]
  tau <- rep(NA, length(tau_str))
  for(i in 1:length(tau_str)){
    tau[i] <- eval(parse(text = tau_str[i]))
  }
  
  # Generate the zeta vector
  if(length(opt$zetas) == 0){
    zetas <- NULL
  }else{
    zetas <- c(as.double(strsplit(opt$zetas, ',')[[1]]))
  }
  
  # Generate the true x values
  set.seed(opt$mc)
  
  x <- runif(opt$nsamp, min = 0, max = 40)
  X <- calc_design_matrix(x, opt$power, zetas)
  
  # Calculate the measurement error variance as a function of x
  sigma2_q_l <- opt$delta_l*x^2
  sigma2_q_u <- opt$delta_u*x^2
  sigma2_q <- (opt$var_weight*sigma2_q_l + (1-opt$var_weight)*sigma2_q_u)
  
  true_signal <- c(X %*% beta)
  sigma2_w_l <- opt$delta_l*true_signal^2
  sigma2_w_u <- opt$delta_u*true_signal^2
  sigma2_w <- (opt$var_weight*sigma2_w_l + (1-opt$var_weight)*sigma2_w_u)
  
  sigma2 <- exp(tau[1] + tau[2]*x)
  
  
  # Generate errors in the model and the measurement errors
  if(opt$dist == 'norm'){
    
    v <- rnorm(opt$nsamp, mean = 0, sd = sqrt(sigma2))
    w <- rnorm(opt$nsamp, mean = 0, sd = sqrt(sigma2_w))
    q <- rnorm(opt$nsamp, mean = 0, sd = sqrt(sigma2_q))
    
    
  } else if(opt$dist == 't'){
    
    t_tc_df <- opt$dist_param
    
    chi_rv <- rchisq(n = opt$nsamp, df = t_tc_df)
    
    norm_rv1 <- rnorm(opt$nsamp, mean = 0, sd = sqrt(((t_tc_df-2)/t_tc_df)*sigma2))
    norm_rv2 <- rnorm(opt$nsamp, mean = 0, sd = sqrt(((t_tc_df-2)/t_tc_df)*sigma2_w))
    norm_rv3 <- rnorm(opt$nsamp, mean = 0, sd = sqrt(((t_tc_df-2)/t_tc_df)*sigma2_q))
    
    v <- norm_rv1/sqrt(chi_rv/t_tc_df)
    w <- norm_rv2/sqrt(chi_rv/t_tc_df)
    q <- norm_rv3/sqrt(chi_rv/t_tc_df)
    
  } else if(opt$dist == 'gnorm'){
    # From Gomez 2007 propositions 3.1 and 3.2 (see also Pascal 2013)
    
    v <- rep(NA, opt$nsamp)
    w <- rep(NA, opt$nsamp)
    q <- rep(NA, opt$nsamp)
    
    gn_tc_theta <- opt$dist_param
    
    for(i in 1:opt$nsamp){
      
      mu <- rep(0, 3)
      var_scale <- ((2^(1/gn_tc_theta))*gamma(5/(2*gn_tc_theta)))/(3*gamma(3/(2*gn_tc_theta)))
      Sigma <- diag(c(sigma2[i], sigma2_w[i], sigma2_q[i]))/var_scale
      tau <- rgamma(1, shape = 3/(2*gn_tc_theta), scale = 2)^(1/(2*gn_tc_theta))
      U <- rnorm(3)
      U_normalized <- U/sqrt(sum(U^2))
      
      mggd_sample <- c(tau * sqrt(Sigma) %*% U_normalized)
      #note: can only use sqrt(Sigma) this way because it's diagonal and p.d.
      
      v[i] <- mggd_sample[1]
      w[i] <- mggd_sample[2]
      q[i] <- mggd_sample[3]
    }
    
  }
  
  
  # Generate the observed x and y values
  truey <- c(X %*% beta + v)
  y_obs <- truey + w
  x_obs <- x + q
  
  # Generate the observed design matrix
  Psi_Xobs <- calc_design_matrix(x_obs, opt$power, zetas)
  
  
  # Split into observed (training) data and new data that we need to calculate prediction intervals for
  train_idx <- sample(1:opt$nsamp, opt$nobs)
  test_idx <- setdiff(1:opt$nsamp, train_idx)
  
  sigma2_q_l_train <- sigma2_q_l[train_idx]
  sigma2_q_u_train <- sigma2_q_u[train_idx]
  sigma2_w_train <- sigma2_w[train_idx]
  sigma2_q_train <- sigma2_q[train_idx]
  y_obs_train <- y_obs[train_idx]
  Psi_Xobs_train <- Psi_Xobs[train_idx,]
  truex_train <- x[train_idx]
  truey_train <- truey[train_idx]
  
  sigma2_q_l_test <- sigma2_q_l[test_idx]
  sigma2_q_u_test <- sigma2_q_u[test_idx]
  sigma2_w_test <- sigma2_w[test_idx]
  sigma2_q_test <- sigma2_q[test_idx]
  y_obs_test <- y_obs[test_idx]
  Psi_Xobs_test <- Psi_Xobs[test_idx,]
  truex_test <- x[test_idx]
  truey_test <- truey[test_idx]
  
  
  # Return the data that we need for the simulations
  return_vals <- list(sigma2_q_l_train, sigma2_q_u_train,
                      sigma2_q_l_train, sigma2_q_u_train,
                      Psi_Xobs_train, y_obs_train,
                      
                      sigma2_q_l_test, sigma2_q_u_test,
                      sigma2_q_l_test, sigma2_q_u_test,
                      Psi_Xobs_test, y_obs_test,
                      
                      truex_train, truex_test,
                      truey_train, truey_test,
                      zetas, beta, sigma2)
  
  names(return_vals) <- c('sigma2_w_l_train', 'sigma2_w_u_train',
                          'sigma2_q_l_train', 'sigma2_q_u_train',
                          'Psi_Xobs_train', 'y_obs_train',
                          
                          'sigma2_w_l_test', 'sigma2_w_u_test',
                          'sigma2_q_l_test', 'sigma2_q_u_test',
                          'Psi_Xobs_test', 'y_obs_test',
                          
                          'truex_train', 'truex_test',
                          'truey_train', 'truey_test',
                          'zetas', 'true_beta', 'sigma2')
  
  return(return_vals)
  
}








####################################
# Functions for estimating the regression parameters
####################################


calc_spline_part <- function(x_i, zeta_j){
  if(x_i > zeta_j){
    spline_part <- x_i - zeta_j
  }else{
    spline_part <- 0
  }
  return(spline_part)
}

calc_nablaPsi_xi <- function(x_i, zeta, p){
  if(is.null(zeta)){
    K <- 0
  }else{
    K <- length(zeta)
  }
  m <- p + K
  
  nablaPsi_xi <- rep(NA, m)
  for(j in 1:p){
    nablaPsi_xi[j] <- j*x_i^(j-1)
  }
  
  if(K > 0){
    for(j in (p+1):m){
      k <- j - p
      nablaPsi_xi[j] <- + p*(calc_spline_part(x_i, zeta[k])^(p-1))*1*(x_i > zeta[k])
    }
  }
  
  return(nablaPsi_xi)
}



calc_beta_hat_mc <- function(x, Psi_Xobs, y_obs, zetas, sigma_q, power){
  
  # Calculate Sigma_q
  n <- length(sigma_q)
  mean_sigmaq <- diag(length(zetas)+1)*0
  for(i in 1:n){
    nablaPsi_xi <- calc_nablaPsi_xi(x[i], zetas, power)
    
    mean_sigmaq <- mean_sigmaq + outer(nablaPsi_xi, nablaPsi_xi)*sigma_q[i]
    
  }
  
  vector0 <- matrix(0, nrow = length(zetas) + 1, ncol = 1)
  mean_Sigmaq <- rbind(cbind(0, t(vector0)), cbind(vector0, mean_sigmaq))
  
  # Calculate the moment corrected MC
  denom <- solve((t(Psi_Xobs) %*% Psi_Xobs) - mean_Sigmaq)
  num <- t(Psi_Xobs) %*% y_obs
  beta_hat <- denom %*% num
  
  return(beta_hat)
}









####################################
# Functions for estimating the true covariate values
####################################

calc_deltaS__delta_xi <- function(testx, p, zeta,
                                  y_i_obs, x_i_obs, beta_hat,
                                  sigma2_p_sigma2w_i, sigma2_q_i){
  
  Psi_xi <- calc_design_matrix(testx, p, zeta)
  nabla_Psi_xi <- c(0, calc_nablaPsi_xi(testx, zeta, p))
  y_piece <- as.numeric((2*(y_i_obs - c(Psi_xi %*% beta_hat))*(- nabla_Psi_xi %*% beta_hat))/sigma2_p_sigma2w_i)
  x_piece <- 2*(x_i_obs - testx)/sigma2_q_i
  
  return(y_piece + x_piece)
  
}

# NOTE: ASSUMES THAT THE POWER IS 1
calc_xhat_i <- function(y_i_obs, sigma2_p_sigma2w_i,
                        Psi_Xobs_idot, sigma2_q_i,
                        beta_hat, zeta, epsilon = 0.01){
  
  # Initialize the vector of potential x_hats
  x_hats <- c()
  
  # Set x_i_obs
  x_i_obs <- Psi_Xobs_idot[2]
  
  if(sum(is.null(zeta))>0){
    y_res <- (y_i_obs - beta_hat[1] - beta_hat[2]*x_i_obs)
    num <- beta_hat[2]*sigma2_q_i*y_res
    denom <- sigma2_p_sigma2w_i + sigma2_q_i*beta_hat[2]^2
    
    x_hat_final <- x_i_obs + num/denom
  } else{
    
    
    # Which region is the true x in?
    # look at derivates of the loss function
    # Loss function is (y_obs - Psi(X)'beta)'(y - Psi(X)'beta) + (x_obs - x)'(x_obs - x)
    # Need derivative to go from negative at start of region to positive at end of region
    
    # For the first region, we just need to check that the derivative is positive
    # at the first knot point
    deltaS__delta_xi <- calc_deltaS__delta_xi(zeta[1], 1, zeta,
                                              y_i_obs, x_i_obs, beta_hat,
                                              sigma2_p_sigma2w_i, sigma2_q_i)
    
    if(deltaS__delta_xi >= 0){
      
      y_res <- (y_i_obs - beta_hat[1] - beta_hat[2]*x_i_obs)
      num <- beta_hat[2]*sigma2_q_i*y_res
      denom <- sigma2_p_sigma2w_i + sigma2_q_i*beta_hat[2]^2
      
      x_hat <- x_i_obs + num/denom
      x_hats <- append(x_hats, x_hat)
      
    }
    
    # For the last region, we need to check
    # the knot point (negative at knot and positive just after)
    deltaS__delta_xi_left <- calc_deltaS__delta_xi(zeta[length(zeta)], 1, zeta,
                                                   y_i_obs, x_i_obs, beta_hat,
                                                   sigma2_p_sigma2w_i, sigma2_q_i)
    deltaS__delta_xi_right <- calc_deltaS__delta_xi(zeta[length(zeta)] + epsilon, 1, zeta,
                                                    y_i_obs, x_i_obs, beta_hat,
                                                    sigma2_p_sigma2w_i, sigma2_q_i)
    
    if(deltaS__delta_xi_left <= 0 & deltaS__delta_xi_right >= 0){
      x_hats <- append(x_hats, zeta[length(zeta)])
    }
    # or that the derivative is negative at the last knot point
    deltaS__delta_xi <- calc_deltaS__delta_xi(zeta[length(zeta)], 1, zeta,
                                              y_i_obs, x_i_obs, beta_hat,
                                              sigma2_p_sigma2w_i, sigma2_q_i)
    
    if(deltaS__delta_xi <= 0){
      
      y_res <- c(y_i_obs - Psi_Xobs_idot %*% beta_hat)
      nabla_Psi_zetaep <- calc_nablaPsi_xi(zeta[length(zeta)] + epsilon, zeta, 1)
      sum_betas <- as.numeric(c(0, nabla_Psi_zetaep) %*% beta_hat)
      num <-  sum_betas*sigma2_q_i*y_res
      denom <- sigma2_p_sigma2w_i + sigma2_q_i*sum_betas^2
      
      x_hat <- x_i_obs + num/denom
      x_hats <- append(x_hats, x_hat)
      
    }
    
    # For regions in between
    for(j in 1:(length(zeta) - 1)){
      # Check the knot point
      deltaS__delta_xi_left <- calc_deltaS__delta_xi(zeta[j], 1, zeta,
                                                     y_i_obs, x_i_obs, beta_hat,
                                                     sigma2_p_sigma2w_i, sigma2_q_i)
      deltaS__delta_xi_right <- calc_deltaS__delta_xi(zeta[j] + epsilon, 1, zeta,
                                                      y_i_obs, x_i_obs, beta_hat,
                                                      sigma2_p_sigma2w_i, sigma2_q_i)
      
      if(deltaS__delta_xi_left <= 0 & deltaS__delta_xi_right >= 0){
        x_hats <- append(x_hats, zeta[j])
      }
      
      # Check the region
      deltaS__delta_xi_left <- calc_deltaS__delta_xi(zeta[j] + epsilon, 1, zeta,
                                                     y_i_obs, x_i_obs, beta_hat,
                                                     sigma2_p_sigma2w_i, sigma2_q_i)
      deltaS__delta_xi_right <- calc_deltaS__delta_xi(zeta[j + 1], 1, zeta,
                                                      y_i_obs, x_i_obs, beta_hat,
                                                      sigma2_p_sigma2w_i, sigma2_q_i)
      
      if(deltaS__delta_xi_left <= 0 & deltaS__delta_xi_right >= 0){
        y_res <- c(y_i_obs - Psi_Xobs_idot %*% beta_hat)
        nabla_Psi_zetaep <- calc_nablaPsi_xi(zeta[j] + epsilon, zeta, 1)
        sum_betas <- as.numeric(c(0, nabla_Psi_zetaep) %*% beta_hat)
        num <-  sum_betas*sigma2_q_i*y_res
        denom <- sigma2_p_sigma2w_i + sigma2_q_i*sum_betas^2
        
        x_hat <- x_i_obs + num/denom
        x_hats <- append(x_hats, x_hat)
      }
    }
    
    
    
    # If there are multiple x-hats, choose the one that minimizes the loss function
    eligible_xhats <- x_hats
    # if none of the zeros are eligible set the x-hat to the observed
    if(length(eligible_xhats) == 0){
      eligible_xhats <- x_i_obs
    }
    nx_hats <- length(eligible_xhats)
    
    
    if(nx_hats > 1){
      S_xis <- rep(NA, nx_hats)
      for(j in 1:nx_hats){
        Psi_xi_hat <- calc_design_matrix(eligible_xhats[j], 1, zeta)
        mu_xi <- c(Psi_xi_hat %*% beta_hat)
        S_xis[j] <- ((x_i_obs - eligible_xhats[j])^2)/sigma2_q_i + 
          ((y_i_obs - mu_xi)^2)/(sigma2_p_sigma2w_i)
      }
      x_hat_final <- eligible_xhats[which(S_xis == min(S_xis))]
    }else{
      x_hat_final <- eligible_xhats
    }
  }
  
  return(x_hat_final)
}



####################################
# Function that estimates beta and eta together
####################################


calc_beta_eta_hat <- function(Psi_Xobs, y_obs, sigma2_w, sigma2_q,
                              sigma2_w_l, sigma2_w_u, zeta,
                              tol = 0.01, max_iter = 100){
  
  
  # Find the number of observations
  n <- length(y_obs)
  
  # Step 1: Initialize estimate of true x value
  xhats <- Psi_Xobs[,2]
  
  # Step 2: Initialize beta_hat using the moment-corrected estimator
  beta_hat <- calc_beta_hat_mc(xhats, Psi_Xobs, y_obs,
                               zeta, sigma2_q, 1)
  yhat <- c(Psi_Xobs %*% beta_hat)
  
  # Step 3: Initialize intrinsic sigma2_plus_sigma2_w_hat
  sigma2_plus_sigma2_w_hat <- sigma2_w + 1
  
  
  # Repeat steps until convergence
  diff <- 3*tol
  iter <- 0
  
  while((diff > tol) & (iter <= max_iter)){
    
    # Step 4: Estimate the true x's using least squares
    # use the current estimate of beta, sigma2_hat + sigma2_w, and sigma2_q
    xhats <- rep(NA, n)
    for(i in 1:n){
      xhats[i] <- calc_xhat_i(y_obs[i], sigma2_plus_sigma2_w_hat[i],
                              Psi_Xobs[i,], sigma2_q[i],
                              beta_hat, zeta)
    }
    
    
    # Step 5: Estimate the variances for the observed ys
      # enforce lower bound for sigma2_plus_sigma2_w_hat
    Psi_Xhats <- calc_design_matrix(xhats, 1, zeta)
    y_hat <- c(Psi_Xhats %*% beta_hat)
    sigma2_plus_sigma2_w_hat <- (y_obs - y_hat)^2
    too_small <- (sigma2_plus_sigma2_w_hat < sigma2_w_l)
    sigma2_plus_sigma2_w_hat[too_small] <- sigma2_w_l[too_small]
    
    
    # Step 6: Estimate beta using estimate of eta
    eta_hat <- rep(NA, n)
    for(i in 1:n){
      nablaPsi_xi_hat <- calc_nablaPsi_xi(xhats[i], zeta, 1)
      beta_nablaPsi_xi <- c(nablaPsi_xi_hat %*% beta_hat[2:length(beta_hat)])
      
      eta_hat[i] <- sigma2_plus_sigma2_w_hat[i] + sigma2_q[i]*beta_nablaPsi_xi^2
    }
    
    Sigma_eta_inv <- diag(1/eta_hat)
    
    Xt_Sigmainv_X__inv <- solve(t(Psi_Xobs) %*% Sigma_eta_inv %*% Psi_Xobs)
    Xt_Sigmainv_y <- t(Psi_Xobs) %*% Sigma_eta_inv %*% y_obs
    
    beta_hat_mp1 <- Xt_Sigmainv_X__inv %*% Xt_Sigmainv_y
    yhat_mp1 <- c(Psi_Xobs %*% beta_hat_mp1)
    
    # Update iteration parameters
    diff <- sum((yhat_mp1 - yhat)^2)/sum((yhat)^2)
    
    beta_hat <- beta_hat_mp1
    yhat <- yhat_mp1
    iter <- iter + 1
  }
  
  if((iter-1) == max_iter){
    warning('Max iter reached')
  }
  
  return_list <- list(beta_hat, iter - 1, diff)
  names(return_list) <- c('beta_hat', 'iter', 'diff')
  
  return(return_list)
  
}


####################################
# Function that calculates the non-conformity score for one removed observation
####################################

calc_ncss <- function(rm_idx, Psi_Xobs_1test, y_obs_1test,
                      sigma2_w_hat_1test, sigma2_q_hat_1test,
                      sigma2_w_l_1test, sigma2_w_u_1test, zeta,
                      beta_tol = 0.01, beta_max_iter = 100,
                      xtol = 0.01, xmax_iter = 100){
  
  
  # a) Remove the observation
  mi_Psi_Xobs <- Psi_Xobs_1test[-rm_idx,]
  mi_y_obs <- y_obs_1test[-rm_idx]
  mi_sigma2_w <- sigma2_w_hat_1test[-rm_idx]
  mi_sigma2_q <- sigma2_q_hat_1test[-rm_idx]
  mi_sigma2_w_l <- sigma2_w_l_1test[-rm_idx]
  mi_sigma2_w_u <- sigma2_w_u_1test[-rm_idx]
  
  # b) Calculate beta_hat without that observation
  calc_beta_eta_hat_res <- calc_beta_eta_hat(mi_Psi_Xobs, mi_y_obs, mi_sigma2_w, mi_sigma2_q,
                                             mi_sigma2_w_l, mi_sigma2_w_u, zeta)
  beta_hat_mi <- calc_beta_eta_hat_res$beta_hat
  
  # c) Estimate the \hat{x}_i using the final estimate of beta
  xtol = 0.01
  xmax_iter = 100
  # Initialize xhat
  xhat <- calc_xhat_i(y_obs_1test[rm_idx], sigma2_w_hat_1test[rm_idx],
                      Psi_Xobs_1test[rm_idx,], sigma2_q_hat_1test[rm_idx],
                      beta_hat_mi, zeta)
  
  
  xdiff <- 3*xtol
  xiter <- 0
  while((xdiff > xtol) & (xiter < xmax_iter)){
    
    # Update sigma2_plus_sigma2_w_hat; always in bounds
    Psi_xi_hat <- calc_design_matrix(xhat, 1, zeta)
    sigma2_plus_sigma2_w_hat <- c((y_obs_1test[rm_idx] - Psi_xi_hat %*% beta_hat_mi)^2)
    # enforce lower bound
    if(sigma2_plus_sigma2_w_hat < sigma2_w_l_1test[rm_idx]){
      sigma2_plus_sigma2_w_hat <- sigma2_w_l_1test[rm_idx]
    }
    
    # Update x_hat
    xhat_mp1 <- calc_xhat_i(y_obs_1test[rm_idx], sigma2_plus_sigma2_w_hat,
                            Psi_Xobs_1test[rm_idx,], sigma2_q_hat_1test[rm_idx],
                            beta_hat_mi, zeta)
    
    # Update diff, iter, and xhat
    xdiff <- (xhat_mp1 - xhat)^2
    xiter <- xiter + 1
    xhat <- xhat_mp1
  }
  
  
  # Use the final values of x_i_hat, sigma2_plus_sigma2_w_hat, and sigma2_q_hat
  # to calculate the non-conformity score of the observation
  nabla_Psi_xi_hat <- calc_nablaPsi_xi(xhat, zeta, 1)
  beta_nablaPsi_xi <- c(nabla_Psi_xi_hat %*% beta_hat_mi[2:length(beta_hat_mi)])
  
  eta_hat <- c(sigma2_plus_sigma2_w_hat + sigma2_q_hat_1test[rm_idx]*beta_nablaPsi_xi^2)
  num <- c(y_obs_1test[rm_idx] - Psi_Xobs_1test[rm_idx,] %*% beta_hat_mi)^2
  r_idx_obs <- num/eta_hat
  
  return_vals <- list(r_idx_obs, xhat)
  names(return_vals) <- c('r_idx_obs', 'xhat')
  
  return(return_vals)
}