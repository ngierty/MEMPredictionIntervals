################################################################################################
# Libraries and functions
################################################################################################

library(optparse)
source('functions.R')
library(ggplot2)
library(latex2exp)

# Extract the parameters from the bash file
option_list = list(
  make_option('--mc', type='integer', default=1,
              help='Test Index')
);

opt_params = OptionParser(option_list=option_list)
opt = parse_args(opt_params)

###############################################
# Calculate the non-conformity scores
###############################################

# Load in the data for calculating the non-conformity scores
load(file = paste0('inter/mr_data/Psi_Xobs_1test__test_idx=', opt$mc, '.RData'))
load(file = paste0('inter/mr_data/y_obs_1test__test_idx=', opt$mc, '.RData'))
load(file = paste0('inter/mr_data/sigma2_w_l_1test__test_idx=', opt$mc, '.RData'))
load(file = paste0('inter/mr_data/sigma2_w_u_1test__test_idx=', opt$mc, '.RData'))
load(file = paste0('inter/mr_data/sigma2_w_hat_1test__test_idx=', opt$mc, '.RData'))
load(file = paste0('inter/mr_data/sigma2_q_hat_1test__test_idx=', opt$mc, '.RData'))

# Set the break points based on Ma (2019) and Ma (2021)
zetas <- c(1.5, 3.84, 8, 13)
log10_zetas <- log10(zetas)


# Calculate the non-conformity scores of all the observations (plus the new test one)
calc_ncss_vals <- sapply(1:nrow(Psi_Xobs_1test), calc_ncss,
                         Psi_Xobs_1test, y_obs_1test,
                         sigma2_w_hat_1test, sigma2_q_hat_1test,
                         sigma2_w_l_1test, sigma2_w_u_1test, log10_zetas)

save(calc_ncss_vals, file = paste0('inter/mr_data/calc_ncss_vals__test_idx=', opt$mc, '.RData'))



























