####################################################################################
# GOAL: Run one simulation where we determine whether the conformal prediction 
#       covers the true value of the response
#           Note, we'll use different test points for the different nominal 
#             coverages we're interested in
#           We'll combine multiple simulations to determine the empirical coverage probability
####################################################################################

################################################
# Load libraries, functions, and extract parameters from bash file
################################################

library(optparse)
library(ggplot2)
library(latex2exp)
source('functions.R')

# Extract the parameters from the bash file
option_list = list(
  make_option('--mc', type='integer', default=NULL,
              help='Monte Carlo Seed'),
  make_option('--tau', type='character', default='log(0.5), 0',
              help='Variance of the intrinsic scatter in equation (1) of the manuscript [default %default]'),
  make_option('--beta', type='character', default='1,3',
              help='Regression parameter vector as string with params separated by comma [default c(%default)]'),
  make_option('--power', type='integer', default=1, help='Power of the splines'),
  make_option('--zetas', type='character', default=NULL, 
              help='Values of the spline knots'),
  make_option('--delta_l', type = 'double', default = 0.06,
              help='Lower bound of variance of measurement error as percentage of x [default %default]'),
  make_option('--delta_u', type = 'double', default = 0.07,
              help='Upper bound of variance of measurement error as percentage of x [default %default]'),
  make_option('--var_weight', type = 'double', default = 0.3,
              help='Weight of variance bounds'),
  make_option('--nobs', type='integer', default=100, help='Sample size'),
  make_option('--dist', type='character', default='norm',
              help='Measurement error distribution; options: norm, gnorm, cauchy, t'),
  make_option('--dist_param', type='double', default=NA,
              help='Degrees of freedom for t or theta for generalized normal')
);

opt_params = OptionParser(option_list=option_list)
opt = parse_args(opt_params)

type1_err <- seq(0.05, 0.95, by = 0.05)
nerr <- length(type1_err)

opt$nsamp <- opt$nobs + nerr


################################################
# Load libraries, functions, and extract parameters from bash file
################################################

# Generate the data
me_data <- gen_me_data(opt)


var_names <- names(me_data)

# Extract the data from the list
for(i in 1:length(var_names)){
  eval(parse(text = paste(var_names[i], ' <- me_data$', var_names[i], sep = '')))  
}


# Save the truex values for calculating the empirical coverages
save(truex_test, file = paste0('inter/sim_data/truex_test__mc=', opt$mc,
                               '__beta=', opt$beta, '__tau=', opt$tau,
                               '__zetas=', opt$zetas, '__nobs=', opt$nobs,
                               '__dist=', opt$dist, '.RData'))

# Calculate the midpoints of the known bounds, same way that astronomers would
sigma2_w_hat_train <- (sigma2_w_l_train + sigma2_w_u_train)/2
sigma2_q_hat_train <- (sigma2_q_l_train + sigma2_q_u_train)/2

sigma2_w_hat_test <- (sigma2_w_l_test + sigma2_w_u_test)/2
sigma2_q_hat_test <- (sigma2_q_l_test + sigma2_q_u_test)/2

################################################
# Calculate the non-conformity scores
  # these are the deleted observed residuals

    # For each new test value
    # 1) Append the new test value to the training data
          # Then for each observation
          # a) Remove the observation
          # b) Calculate beta_hat without that observation
          # c) \hat{x} using the final estimate of beta
          # d) Calculate the non-conformity score, equation (2) of the manuscript,
                # using beta_hat(-i) and eta_hat(-i)
    # 2) Save the non-conformity scores
################################################

st <- Sys.time()
ncs_pvalues <- rep(NA, nerr)
test_xhats <- rep(NA, nerr)
ncs_region_pvalue <- rep(NA, nerr)

for(test_idx in 1:nerr){

  # Append the new test value to the training data
  Psi_Xobs_1test <- rbind(Psi_Xobs_train, Psi_Xobs_test[test_idx,])
  y_obs_1test <- append(y_obs_train, y_obs_test[test_idx])
  sigma2_w_hat_1test <- append(sigma2_w_hat_train, sigma2_w_hat_test[test_idx])
  sigma2_w_l_1test <- append(sigma2_w_l_train, sigma2_w_l_test[test_idx])
  sigma2_w_u_1test <- append(sigma2_w_u_train, sigma2_w_u_test[test_idx])
  sigma2_q_hat_1test <- append(sigma2_q_hat_train, sigma2_q_hat_test[test_idx])
  
  
  # Calculate the non-conformity scores of all the observations (plus the new test one)
  calc_ncss_vals <- sapply(1:(opt$nobs + 1), calc_ncss, Psi_Xobs_1test, y_obs_1test,
                           sigma2_w_hat_1test, sigma2_q_hat_1test,
                           sigma2_w_l_1test, sigma2_w_u_1test, zetas)
  
  # Calculate the p-value of the non-conformity score
  ncss <- unlist(calc_ncss_vals[1,])
  ncs_pvalues[test_idx] <- mean(ncss > ncss[opt$nobs + 1])
  
  
  # Calculate the p-value of the non-conformity score based on the estimated region
  xhats <- unlist(calc_ncss_vals[2,])
    # Save the xhat value of the test point
  test_xhats[test_idx] <- xhats[opt$nobs + 1]
  
  estx_region1_flag <- (xhats < zetas[1])
  estx_region2_flag <- (zetas[1] <= xhats) & (xhats < zetas[2])
  estx_region3_flag <- (xhats >= zetas[2])
  
  if(estx_region1_flag[opt$nobs + 1]){
    ncs_region_pvalue[test_idx] <- mean(ncss[estx_region1_flag] > ncss[opt$nobs + 1])
  }else if(estx_region2_flag[opt$nobs + 1]){
    ncs_region_pvalue[test_idx] <- mean(ncss[estx_region2_flag] > ncss[opt$nobs + 1])
  }else{
    ncs_region_pvalue[test_idx] <- mean(ncss[estx_region3_flag] > ncss[opt$nobs + 1])
  }
  
  
}

warnings()

save(ncs_pvalues, file = paste0('inter/sim_data/ncs_pvalues__mc=', opt$mc,
                                '__beta=', opt$beta, '__tau=', opt$tau,
                                '__zetas=', opt$zetas, '__nobs=', opt$nobs,
                                '__dist=', opt$dist, '.RData'))


save(test_xhats, file = paste0('inter/sim_data/test_xhats__mc=', opt$mc,
                          '__beta=', opt$beta, '__tau=', opt$tau,
                          '__zetas=', opt$zetas, '__nobs=', opt$nobs,
                          '__dist=', opt$dist, '.RData'))


save(ncs_region_pvalue, file = paste0('inter/sim_data/ncs_region_pvalue__mc=', opt$mc,
                                '__beta=', opt$beta, '__tau=', opt$tau,
                                '__zetas=', opt$zetas, '__nobs=', opt$nobs,
                                '__dist=', opt$dist, '.RData'))


rt <- Sys.time() - st
print(rt)






























