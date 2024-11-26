################################################################################################
# Setup data to perform non-conformal prediction
################################################################################################

library(ggplot2)
library(latex2exp)
source('functions.R')

########################################
# Clean up the data
########################################

# Read in the data
mr <- read.csv('raw/NASA_exoplanet_archive.csv', skip = 21, header = TRUE)

mr_bmm <- read.csv('raw/NASA_MRF_confirmed.csv', header = TRUE)

# Remove observations that we don't care about
  # Remove observations that are likely brown dwarves and other stars
planet <- !is.na(mr$pl_bmassj) & (mr$pl_bmassj < 13)
pos_rad <- !is.na(mr$pl_rade)
pos_mass <- !is.na(mr$pl_bmasse)
default_param <- (mr$default_flag == 1) # keeps only those that are the default planetary parameters
not_contr <- (mr$pl_controv_flag == 0) # removes disputed planets

# clean up the error standard deviations
n <- nrow(mr)
mr$sigma2_m_l <- rep(NA, n)
mr$sigma2_m_u <- rep(NA, n)
mr$sigma2_r_l <- rep(NA, n)
mr$sigma2_r_u <- rep(NA, n)

for(i in 1:n){
  mr$sigma2_m_l[i] <- min(mr$pl_bmasseerr1[i]^2, abs(mr$pl_bmasseerr2[i])^2)
  mr$sigma2_m_u[i] <- max(mr$pl_bmasseerr1[i]^2, abs(mr$pl_bmasseerr2[i])^2)
  mr$sigma2_r_l[i] <- min(mr$pl_radeerr1[i]^2, abs(mr$pl_radeerr2[i])^2)
  mr$sigma2_r_u[i] <- max(mr$pl_radeerr1[i]^2, abs(mr$pl_radeerr2[i])^2)
}

# use the midpoints as the known values
mr$sigma2_m_hat <- (mr$sigma2_m_l + mr$sigma2_m_u)/2
mr$sigma2_r_hat <- (mr$sigma2_r_l + mr$sigma2_r_u)/2

rad_snr <- mr$pl_rade/sqrt(mr$sigma2_r_hat)
mass_snr <- mr$pl_bmasse/sqrt(mr$sigma2_m_hat)

rad_snr_high <- !is.na(rad_snr) & (rad_snr >= 3)
mass_snr_high <- !is.na(mass_snr) & (mass_snr >= 3)

# keep only the data that we need (based on 'Model sampling.ipynb' from MA 2021)
in_ma2021 <- mr$pl_name %in% mr_bmm$pl_name
mr_keep <- mr[planet & pos_rad & pos_mass & default_param & not_contr
              & rad_snr_high & mass_snr_high, ]

# did any planets get dropped
in_ma2021 <- unique(mr_keep$pl_name) %in% mr_bmm$pl_name
sum(in_ma2021)
    

  # log base 10 transform everything
mr_keep$l10_R <- log10(mr_keep$pl_rade)
mr_keep$l10_M <- log10(mr_keep$pl_bmasse)
mr_keep$l10_sigma2_m_l <- mr_keep$sigma2_m_l/log(10)
mr_keep$l10_sigma2_m_u <- mr_keep$sigma2_m_u/log(10)
mr_keep$l10_sigma2_r_l <- mr_keep$sigma2_r_l/log(10)
mr_keep$l10_sigma2_r_u <- mr_keep$sigma2_r_u/log(10)
mr_keep$l10_sigma2_m_hat <- mr_keep$sigma2_m_hat/log(10)
mr_keep$l10_sigma2_r_hat <- mr_keep$sigma2_r_hat/log(10)


########################################
# Split the data into train and test sets
########################################

# Subset the data into training and test sets
set.seed(9586)
n <- nrow(mr_keep)
ntest <- 300
ntrain <- n - ntest

train_idx <- sample(1:n, ntrain)
test_idx <- setdiff(1:n, train_idx)


mr_train <- mr_keep[train_idx,]
mr_test <- mr_keep[test_idx,]

# Set the break points based on Ma (2019) and Ma (2021)
zetas <- c(1.5, 3.84, 8, 13)
log10_zetas <- log10(zetas)


# Plot the training data
plot_df <- mr_train
plot_df$l10_sigma_r_hat <- sqrt(plot_df$l10_sigma2_r_hat)
plot_df$l10_sigma_m_hat <- sqrt(plot_df$l10_sigma2_m_hat)
ggplot(plot_df) + geom_point(aes(x = l10_R, y = l10_M)) +
  theme(text = element_text(size = 20)) +
  labs(x = TeX(r'($\log_{10}(R_{\oplus})$)'),
       y = TeX(r'($\log_{10}(M_{\oplus})$)'))
ggsave(paste0('output/mr_data/train_data_scatter.pdf'),
       height = 8, width = 10, units = 'in', dpi = 600)



####################
# Format the data for calculating the non-conformity scores
####################


# Train set
Psi_Xobs_train <- calc_design_matrix(mr_train$l10_R, 1, log10_zetas)
y_obs_train <- mr_train$l10_M
sigma2_w_l_train <- mr_train$l10_sigma2_m_l
sigma2_w_u_train <- mr_train$l10_sigma2_m_u
sigma2_w_hat_train <- mr_train$l10_sigma2_m_hat
sigma2_q_hat_train <- mr_train$l10_sigma2_r_hat


# Test set
Psi_Xobs_test <- calc_design_matrix(mr_test$l10_R, 1, log10_zetas)
y_obs_test <- mr_test$l10_M
sigma2_w_l_test <- mr_test$l10_sigma2_m_l
sigma2_w_u_test <- mr_test$l10_sigma2_m_u
sigma2_w_hat_test <- mr_test$l10_sigma2_m_hat
sigma2_q_hat_test <- mr_test$l10_sigma2_r_hat



#####################
# Output all of the test sets needed to calculate the non-conformity scores
    # the actual non-conformity scores will be computed in parallel on the cluster
        # to reduce computational time
#####################

for(test_idx in 1:ntest){
  
  # Append the new test value to the training data
  Psi_Xobs_1test <- rbind(Psi_Xobs_train, Psi_Xobs_test[test_idx,])
  y_obs_1test <- append(y_obs_train, y_obs_test[test_idx])
  sigma2_w_l_1test <- append(sigma2_w_l_train, sigma2_w_l_test[test_idx])
  sigma2_w_u_1test <- append(sigma2_w_u_train, sigma2_w_u_test[test_idx])
  sigma2_w_hat_1test <- append(sigma2_w_hat_train, sigma2_w_hat_test[test_idx])
  sigma2_q_hat_1test <- append(sigma2_q_hat_train, sigma2_q_hat_test[test_idx])
  
  
  save(Psi_Xobs_1test, file = paste0('inter/mr_data/Psi_Xobs_1test__test_idx=', test_idx, '.RData'))
  save(y_obs_1test, file = paste0('inter/mr_data/y_obs_1test__test_idx=', test_idx, '.RData'))
  save(sigma2_w_l_1test, file = paste0('inter/mr_data/sigma2_w_l_1test__test_idx=', test_idx, '.RData'))
  save(sigma2_w_u_1test, file = paste0('inter/mr_data/sigma2_w_u_1test__test_idx=', test_idx, '.RData'))
  save(sigma2_w_hat_1test, file = paste0('inter/mr_data/sigma2_w_hat_1test__test_idx=', test_idx, '.RData'))
  save(sigma2_q_hat_1test, file = paste0('inter/mr_data/sigma2_q_hat_1test__test_idx=', test_idx, '.RData'))

}
