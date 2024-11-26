############################################################################
# Determine the coverage of the conformal prediction intervals
# We can do this by just calculating the p-value (from Algorithm 1) using the true unobserved response value
# The p-value will tell us whether this response value would have been included in the 
# full conformal prediction interval if we had calculated it
# Save the indicator of whether this value is in the full conformal prediction interval
# In the output.R, we calculate the coverage probability using the average of the indicators across all sims
############################################################################

library(optparse)
source('functions.R')
library(ggplot2)
library(latex2exp)

################################################################
# Parameters
################################################################

type1_err <- seq(0.05, 0.95, by = 0.05)
nom_cov <- 1 - type1_err
nerr <- length(type1_err)
tau_list <- list(c('log(3),0'), c('log(3),0.2'))
dist_list <- c('norm', 'gnorm', 't')

beta <- '1,5,4,7'
zetas_str <- '10,25'
nobs <- 100
n_mcs <- 300

################################################################
# Make the scatter plots of the observed data
################################################################

# Make a representative dataset just with more points for visualization purposes
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

all_obs <- data.frame()

for(dist in dist_list){
  for(tau in tau_list){
    opt$mc <- 35
    opt$beta <- beta
    opt$tau <- tau
    opt$zetas <- zetas_str
    opt$nobs <- 700
    opt$dist <- dist
    opt$nsamp <- opt$nobs + length(type1_err)
    if(dist == 'gnorm'){
      opt$dist_param <- 0.365
    }else if(dist == 't'){
      opt$dist_param <- 3
    }
    
    # Generate the data
    me_data <- gen_me_data(opt)
    var_names <- names(me_data)
    
    # Extract the data from the list
    for(i in 1:length(var_names)){
      eval(parse(text = paste(var_names[i], ' <- me_data$', var_names[i], sep = '')))  
    }
    
    x_obs_train <- Psi_Xobs_train[,2]
    if(tau == 'log(3),0'){
      var_case <- 'Homoscedastic'
    }else{
      var_case <- 'Heteroscedastic'
    }
    df <- data.frame(var_case, dist, truex_train, truey_train, x_obs_train, y_obs_train)
    
    # Save
    all_obs <- rbind(all_obs, df)
  }
}

# Clean up names
all_obs$dist <- factor(all_obs$dist, 
                       labels = c(gnorm = TeX('Gen.\\ Normal'),
                                  norm = 'Normal', t = TeX('$t_3$')))


all_obs$var_case <- factor(all_obs$var_case, 
                           levels = c('Homoscedastic', 'Heteroscedastic'))

# Plot the true x against the true y to visualize the intrinsic scatter
ggplot(all_obs) + geom_point(aes(x = truex_train, y = truey_train)) +
  facet_grid(var_case ~ dist, labeller = label_parsed) +
  theme(text = element_text(size = 20)) +
  labs(x = 'x', y = 'y')
ggsave(paste0('output/sim_data/intrinsic_scatter.pdf'),
       height = 8, width = 16, units = 'in', dpi = 600)


# Plot the true y against the observed y to visualize w
ggplot(all_obs) + geom_point(aes(x = truey_train, y = y_obs_train)) +
  facet_grid(var_case ~ dist, labeller = label_parsed) +
  theme(text = element_text(size = 20)) +
  labs(x = 'y', y = TeX('$y^{obs}$'))
ggsave(paste0('output/sim_data/w.pdf'),
       height = 8, width = 16, units = 'in', dpi = 600)


# Plot the true x against the observed x to visualize q
ggplot(all_obs) + geom_point(aes(x = truex_train, y = x_obs_train)) +
  facet_grid(var_case ~ dist, labeller = label_parsed) +
  theme(text = element_text(size = 20)) +
  labs(x = 'x', y = TeX('$x^{obs}$'))
ggsave(paste0('output/sim_data/q.pdf'),
       height = 8, width = 16, units = 'in', dpi = 600)


# Plot final observed values
ggplot(all_obs) + geom_point(aes(x = x_obs_train, y = y_obs_train)) +
  facet_grid(var_case ~ dist, labeller = label_parsed) +
  theme(text = element_text(size = 20)) +
  labs(x = TeX('$x^{obs}$'), y = TeX('$y^{obs}$'))
ggsave(paste0('output/sim_data/observed_data.pdf'),
       height = 8, width = 16, units = 'in', dpi = 600)



################################################################
# Make the coverage plots
################################################################

# empirical coverage across all regions
all_taus_dists_emp_cov_all_regions <- data.frame()

# empirical coverage across regions (based on true x), p-value un-adjusted in real time
all_taus_dists_emp_cov_truex_region1 <- data.frame()
all_taus_dists_emp_cov_truex_region2 <- data.frame()
all_taus_dists_emp_cov_truex_region3 <- data.frame()

# empirical coverage across regions (based on est x), p-value un-adjusted in real time
all_taus_dists_emp_cov_estx_region1 <- data.frame()
all_taus_dists_emp_cov_estx_region2 <- data.frame()
all_taus_dists_emp_cov_estx_region3 <- data.frame()

# empirical coverage across regions (based on est x), p-value adjusted in real time
all_taus_dists_emp_cov_all_regions_adj <- data.frame()
all_taus_dists_emp_cov_estx_region1_adj <- data.frame()
all_taus_dists_emp_cov_estx_region2_adj <- data.frame()
all_taus_dists_emp_cov_estx_region3_adj <- data.frame()

for(tau in tau_list){
  
  for(dist in dist_list){
    
    #########################################
    # Store all of the non-conformity score p-values
    #########################################
    ncs_pvalues_all <- matrix(NA, nrow = nerr, ncol = n_mcs)
    ncs_pvalues_all_adj <- matrix(NA, nrow = nerr, ncol = n_mcs)
    
    ncs_truex_region1 <- matrix(NA, nrow = nerr, ncol = n_mcs)
    ncs_truex_region2 <- matrix(NA, nrow = nerr, ncol = n_mcs)
    ncs_truex_region3 <- matrix(NA, nrow = nerr, ncol = n_mcs)
    
    ncs_estx_region1 <- matrix(NA, nrow = nerr, ncol = n_mcs)
    ncs_estx_region2 <- matrix(NA, nrow = nerr, ncol = n_mcs)
    ncs_estx_region3 <- matrix(NA, nrow = nerr, ncol = n_mcs)
    
    ncs_estx_region1_adj <- matrix(NA, nrow = nerr, ncol = n_mcs)
    ncs_estx_region2_adj <- matrix(NA, nrow = nerr, ncol = n_mcs)
    ncs_estx_region3_adj <- matrix(NA, nrow = nerr, ncol = n_mcs)
    
    for(mc in 1:n_mcs){
      
      # Get the saved p-values
      ncs_pvalues <- rep(NA, nerr)
      
      load(file = paste0('inter/sim_data/ncs_pvalues__mc=', mc,
                         '__beta=', beta, '__tau=', tau,
                         '__zetas=', zetas_str, '__nobs=', nobs,
                         '__dist=', dist, '.RData'))
      
      # Save all of them ignoring region differences
      ncs_pvalues_all[,mc] <- ncs_pvalues
      
      # Save the p-values by region (determined by true x value)
      zetas_num <- c(as.double(strsplit(zetas_str, ',')[[1]]))
      load(file = paste0('inter/sim_data/truex_test__mc=', mc,
                         '__beta=', beta, '__tau=', tau,
                         '__zetas=', zetas_str, '__nobs=', nobs,
                         '__dist=', dist, '.RData'))
      truex_region1_flag <- (truex_test < zetas_num[1])
      truex_region2_flag <- (zetas_num[1] <= truex_test) & (truex_test < zetas_num[2])
      truex_region3_flag <- (truex_test >= zetas_num[2])
      
      
      ncs_truex_region1[truex_region1_flag, mc] <- ncs_pvalues[truex_region1_flag]
      ncs_truex_region2[truex_region2_flag, mc] <- ncs_pvalues[truex_region2_flag]
      ncs_truex_region3[truex_region3_flag, mc] <- ncs_pvalues[truex_region3_flag]
      
      
      # Save the p-values by region (determined by estimated x value)
      load(file = paste0('inter/sim_data/test_xhats__mc=', mc,
                         '__beta=', beta, '__tau=', tau,
                         '__zetas=', zetas_str, '__nobs=', nobs,
                         '__dist=', dist, '.RData'))
      estx_region1_flag <- (test_xhats < zetas_num[1])
      estx_region2_flag <- (zetas_num[1] <= test_xhats) & (test_xhats < zetas_num[2])
      estx_region3_flag <- (test_xhats >= zetas_num[2])
      
      
      ncs_estx_region1[estx_region1_flag, mc] <- ncs_pvalues[estx_region1_flag]
      ncs_estx_region2[estx_region2_flag, mc] <- ncs_pvalues[estx_region2_flag]
      ncs_estx_region3[estx_region3_flag, mc] <- ncs_pvalues[estx_region3_flag]
      
      
      
      # Save the adjusted p-values by region (determined by estimated x value)
      load(file = paste0('inter/sim_data/ncs_region_pvalue__mc=', mc,
                         '__beta=', beta, '__tau=', tau,
                         '__zetas=', zetas_str, '__nobs=', nobs,
                         '__dist=', dist, '.RData'))
      ncs_pvalues_all_adj[,mc] <- ncs_region_pvalue
      
      ncs_estx_region1_adj[estx_region1_flag, mc] <- ncs_region_pvalue[estx_region1_flag]
      ncs_estx_region2_adj[estx_region2_flag, mc] <- ncs_region_pvalue[estx_region2_flag]
      ncs_estx_region3_adj[estx_region3_flag, mc] <- ncs_region_pvalue[estx_region3_flag]
    }
    
    
    
    #########################################
    # For each type 1 error, calculate the empirical coverage
    #########################################
    
    emp_cov_all_regions <- rep(NA, nerr)
    
    # regions based on true x value
    emp_cov_truex_region1 <- rep(NA, nerr)
    emp_cov_truex_region2 <- rep(NA, nerr)
    emp_cov_truex_region3 <- rep(NA, nerr)
    
    num_truex_region1 <- rep(NA, nerr)
    num_truex_region2 <- rep(NA, nerr)
    num_truex_region3 <- rep(NA, nerr)
    
    # regions based on estimated x value
    emp_cov_estx_region1 <- rep(NA, nerr)
    emp_cov_estx_region2 <- rep(NA, nerr)
    emp_cov_estx_region3 <- rep(NA, nerr)
    
    num_estx_region1 <- rep(NA, nerr)
    num_estx_region2 <- rep(NA, nerr)
    num_estx_region3 <- rep(NA, nerr)
    
    # adjusted p-values by region
    emp_cov_all_regions_adj <- rep(NA, nerr)
    emp_cov_estx_region1_adj <- rep(NA, nerr)
    emp_cov_estx_region2_adj <- rep(NA, nerr)
    emp_cov_estx_region3_adj <- rep(NA, nerr)
    
    
    for(i in 1:nerr){
      
      # All regions
      emp_cov_all_regions[i] <- mean(ncs_pvalues_all[i,] > type1_err[i], na.rm = TRUE)
      emp_cov_all_regions_adj[i] <- mean(ncs_pvalues_all_adj[i,] > type1_err[i], na.rm = TRUE)
      
      # Region 1
      emp_cov_truex_region1[i] <- mean(ncs_truex_region1[i,] > type1_err[i], na.rm = TRUE)
      num_truex_region1[i] <- sum(!is.na(ncs_truex_region1[i,]))
      
      emp_cov_estx_region1[i] <- mean(ncs_estx_region1[i,] > type1_err[i], na.rm = TRUE)
      emp_cov_estx_region1_adj[i] <- mean(ncs_estx_region1_adj[i,] > type1_err[i], na.rm = TRUE)
      num_estx_region1[i] <- sum(!is.na(ncs_estx_region1[i,]))
      
      # Region 2
      emp_cov_truex_region2[i] <- mean(ncs_truex_region2[i,] > type1_err[i], na.rm = TRUE)
      num_truex_region2[i] <- sum(!is.na(ncs_truex_region2[i,]))
      
      emp_cov_estx_region2[i] <- mean(ncs_estx_region2[i,] > type1_err[i], na.rm = TRUE)
      emp_cov_estx_region2_adj[i] <- mean(ncs_estx_region2_adj[i,] > type1_err[i], na.rm = TRUE)
      num_estx_region2[i] <- sum(!is.na(ncs_estx_region2[i,]))
      
      # Region 3
      emp_cov_truex_region3[i] <- mean(ncs_truex_region3[i,] > type1_err[i], na.rm = TRUE)
      num_truex_region3[i] <- sum(!is.na(ncs_truex_region3[i,]))
      
      emp_cov_estx_region3[i] <- mean(ncs_estx_region3[i,] > type1_err[i], na.rm = TRUE)
      emp_cov_estx_region3_adj[i] <- mean(ncs_estx_region3_adj[i,] > type1_err[i], na.rm = TRUE)
      num_estx_region3[i] <- sum(!is.na(ncs_estx_region3[i,]))
    }
    
    #########################################
    # Save the values for plotting purposes later
    #########################################
    
    # all regions
    all_regions_num <- sum(!is.na(ncs_pvalues_all[1,]))
    if(tau == 'log(3),0'){
      var_case <- 'Homoscedastic'
    }else{
      var_case <- 'Heteroscedastic'
    }
    
    # all regions
    df <- data.frame(var_case, dist, nom_cov, emp_cov = emp_cov_all_regions, num = all_regions_num)
    all_taus_dists_emp_cov_all_regions <- rbind(all_taus_dists_emp_cov_all_regions, df)
    
    # region 1, determined by true x
    df <- data.frame(var_case, dist, nom_cov, emp_cov = emp_cov_truex_region1, num = num_truex_region1)
    all_taus_dists_emp_cov_truex_region1 <- rbind(all_taus_dists_emp_cov_truex_region1, df)
    
    # region 2, determined by true x
    df <- data.frame(var_case, dist, nom_cov, emp_cov = emp_cov_truex_region2, num = num_truex_region2)
    all_taus_dists_emp_cov_truex_region2 <- rbind(all_taus_dists_emp_cov_truex_region2, df)
    
    # region 3, determined by true x
    df <- data.frame(var_case, dist, nom_cov, emp_cov = emp_cov_truex_region3, num = num_truex_region3)
    all_taus_dists_emp_cov_truex_region3 <- rbind(all_taus_dists_emp_cov_truex_region3, df)
    
    
    
    
    # region 1, determined by estimated x
    df <- data.frame(var_case, dist, nom_cov, emp_cov = emp_cov_estx_region1, num = num_estx_region1)
    all_taus_dists_emp_cov_estx_region1 <- rbind(all_taus_dists_emp_cov_estx_region1, df)
    
    # region 2, determined by estimated x
    df <- data.frame(var_case, dist, nom_cov, emp_cov = emp_cov_estx_region2, num = num_estx_region2)
    all_taus_dists_emp_cov_estx_region2 <- rbind(all_taus_dists_emp_cov_estx_region2, df)
    
    # region 3, determined by estimated x
    df <- data.frame(var_case, dist, nom_cov, emp_cov = emp_cov_estx_region3, num = num_estx_region3)
    all_taus_dists_emp_cov_estx_region3 <- rbind(all_taus_dists_emp_cov_estx_region3, df)
    
    
    
    
    
    # all regions, adjusted pvalue
    df <- data.frame(var_case, dist, nom_cov, emp_cov = emp_cov_all_regions_adj, num = all_regions_num)
    all_taus_dists_emp_cov_all_regions_adj <- rbind(all_taus_dists_emp_cov_all_regions_adj, df)
    
    # region 1, determined by estimated x
    df <- data.frame(var_case, dist, nom_cov, emp_cov = emp_cov_estx_region1_adj, num = num_estx_region1)
    all_taus_dists_emp_cov_estx_region1_adj <- rbind(all_taus_dists_emp_cov_estx_region1_adj, df)
    
    # region 2, determined by estimated x
    df <- data.frame(var_case, dist, nom_cov, emp_cov = emp_cov_estx_region2_adj, num = num_estx_region2)
    all_taus_dists_emp_cov_estx_region2_adj <- rbind(all_taus_dists_emp_cov_estx_region2_adj, df)
    
    # region 3, determined by estimated x
    df <- data.frame(var_case, dist, nom_cov, emp_cov = emp_cov_estx_region3_adj, num = num_estx_region3)
    all_taus_dists_emp_cov_estx_region3_adj <- rbind(all_taus_dists_emp_cov_estx_region3_adj, df)
  }
}

#########################################
# Add error bounds
#########################################

  # For all regions
all_taus_dists_emp_cov_all_regions$mc_err <- sqrt(nom_cov*(1-nom_cov)/all_taus_dists_emp_cov_all_regions$num)
all_taus_dists_emp_cov_all_regions$lb <- all_taus_dists_emp_cov_all_regions$nom_cov - 
  all_taus_dists_emp_cov_all_regions$mc_err
all_taus_dists_emp_cov_all_regions$ub <- all_taus_dists_emp_cov_all_regions$nom_cov + 
  all_taus_dists_emp_cov_all_regions$mc_err

  # region 1, determined by true x
all_taus_dists_emp_cov_truex_region1$mc_err <- sqrt(nom_cov*(1-nom_cov)/all_taus_dists_emp_cov_truex_region1$num)
all_taus_dists_emp_cov_truex_region1$lb <- all_taus_dists_emp_cov_truex_region1$nom_cov - 
  all_taus_dists_emp_cov_truex_region1$mc_err
all_taus_dists_emp_cov_truex_region1$ub <- all_taus_dists_emp_cov_truex_region1$nom_cov + 
  all_taus_dists_emp_cov_truex_region1$mc_err

  # region 2, determined by true x
all_taus_dists_emp_cov_truex_region2$mc_err <- sqrt(nom_cov*(1-nom_cov)/all_taus_dists_emp_cov_truex_region2$num)
all_taus_dists_emp_cov_truex_region2$lb <- all_taus_dists_emp_cov_truex_region2$nom_cov - 
  all_taus_dists_emp_cov_truex_region2$mc_err
all_taus_dists_emp_cov_truex_region2$ub <- all_taus_dists_emp_cov_truex_region2$nom_cov + 
  all_taus_dists_emp_cov_truex_region2$mc_err

  # region 3, determined by true x
all_taus_dists_emp_cov_truex_region3$mc_err <- sqrt(nom_cov*(1-nom_cov)/all_taus_dists_emp_cov_truex_region3$num)
all_taus_dists_emp_cov_truex_region3$lb <- all_taus_dists_emp_cov_truex_region3$nom_cov - 
  all_taus_dists_emp_cov_truex_region3$mc_err
all_taus_dists_emp_cov_truex_region3$ub <- all_taus_dists_emp_cov_truex_region3$nom_cov + 
  all_taus_dists_emp_cov_truex_region3$mc_err





# region 1, determined by estimated x
all_taus_dists_emp_cov_estx_region1$mc_err <- sqrt(nom_cov*(1-nom_cov)/all_taus_dists_emp_cov_estx_region1$num)
all_taus_dists_emp_cov_estx_region1$lb <- all_taus_dists_emp_cov_estx_region1$nom_cov - 
  all_taus_dists_emp_cov_estx_region1$mc_err
all_taus_dists_emp_cov_estx_region1$ub <- all_taus_dists_emp_cov_estx_region1$nom_cov + 
  all_taus_dists_emp_cov_estx_region1$mc_err

# region 2, determined by estimated x
all_taus_dists_emp_cov_estx_region2$mc_err <- sqrt(nom_cov*(1-nom_cov)/all_taus_dists_emp_cov_estx_region2$num)
all_taus_dists_emp_cov_estx_region2$lb <- all_taus_dists_emp_cov_estx_region2$nom_cov - 
  all_taus_dists_emp_cov_estx_region2$mc_err
all_taus_dists_emp_cov_estx_region2$ub <- all_taus_dists_emp_cov_estx_region2$nom_cov + 
  all_taus_dists_emp_cov_estx_region2$mc_err

# region 3, determined by estimated x
all_taus_dists_emp_cov_estx_region3$mc_err <- sqrt(nom_cov*(1-nom_cov)/all_taus_dists_emp_cov_estx_region3$num)
all_taus_dists_emp_cov_estx_region3$lb <- all_taus_dists_emp_cov_estx_region3$nom_cov - 
  all_taus_dists_emp_cov_estx_region3$mc_err
all_taus_dists_emp_cov_estx_region3$ub <- all_taus_dists_emp_cov_estx_region3$nom_cov + 
  all_taus_dists_emp_cov_estx_region3$mc_err






# adjusted p-values, all regions
all_taus_dists_emp_cov_estx_region3_adj$mc_err <- sqrt(nom_cov*(1-nom_cov)/all_taus_dists_emp_cov_estx_region3_adj$num)
all_taus_dists_emp_cov_estx_region3_adj$lb <- all_taus_dists_emp_cov_estx_region3_adj$nom_cov - 
  all_taus_dists_emp_cov_estx_region3_adj$mc_err
all_taus_dists_emp_cov_estx_region3_adj$ub <- all_taus_dists_emp_cov_estx_region3_adj$nom_cov + 
  all_taus_dists_emp_cov_estx_region3_adj$mc_err


# adjusted p-values region 1, determined by estimated x
all_taus_dists_emp_cov_estx_region1_adj$mc_err <- sqrt(nom_cov*(1-nom_cov)/all_taus_dists_emp_cov_estx_region1_adj$num)
all_taus_dists_emp_cov_estx_region1_adj$lb <- all_taus_dists_emp_cov_estx_region1_adj$nom_cov - 
  all_taus_dists_emp_cov_estx_region1_adj$mc_err
all_taus_dists_emp_cov_estx_region1_adj$ub <- all_taus_dists_emp_cov_estx_region1_adj$nom_cov + 
  all_taus_dists_emp_cov_estx_region1_adj$mc_err

# adjusted p-values region 2, determined by estimated x
all_taus_dists_emp_cov_estx_region2_adj$mc_err <- sqrt(nom_cov*(1-nom_cov)/all_taus_dists_emp_cov_estx_region2_adj$num)
all_taus_dists_emp_cov_estx_region2_adj$lb <- all_taus_dists_emp_cov_estx_region2_adj$nom_cov - 
  all_taus_dists_emp_cov_estx_region2_adj$mc_err
all_taus_dists_emp_cov_estx_region2_adj$ub <- all_taus_dists_emp_cov_estx_region2_adj$nom_cov + 
  all_taus_dists_emp_cov_estx_region2_adj$mc_err

# adjusted p-values region 3, determined by estimated x
all_taus_dists_emp_cov_estx_region3_adj$mc_err <- sqrt(nom_cov*(1-nom_cov)/all_taus_dists_emp_cov_estx_region3_adj$num)
all_taus_dists_emp_cov_estx_region3_adj$lb <- all_taus_dists_emp_cov_estx_region3_adj$nom_cov - 
  all_taus_dists_emp_cov_estx_region3_adj$mc_err
all_taus_dists_emp_cov_estx_region3_adj$ub <- all_taus_dists_emp_cov_estx_region3_adj$nom_cov + 
  all_taus_dists_emp_cov_estx_region3_adj$mc_err









#########################################
# Plot the empirical coverages
#########################################

#####################
# All regions together
#####################

###########
# unadjusted pvalues
###########
all_taus_dists_emp_cov_all_regions$dist <- factor(all_taus_dists_emp_cov_all_regions$dist, 
                                                  labels = c(gnorm = TeX('Gen.\\ Normal'),
                                                             norm = 'Normal', t = TeX('$t_3$')))


all_taus_dists_emp_cov_all_regions$var_case <- factor(all_taus_dists_emp_cov_all_regions$var_case, 
                                                  levels = c('Homoscedastic', 'Heteroscedastic'))

ggplot(all_taus_dists_emp_cov_all_regions) + 
  geom_point(aes(x = nom_cov, y = emp_cov)) + 
  facet_grid(var_case ~ dist, labeller = label_parsed) +
  geom_abline(intercept = 0, slope = 1) + 
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(text = element_text(size = 20)) +
  labs(x = 'Nominal Coverage', y = 'Empirical Coverage')
ggsave(paste0('output/sim_data/all_regions_coverage.pdf'),
       height = 8, width = 16, units = 'in', dpi = 600)




###########
# adjusted pvalues
###########
all_taus_dists_emp_cov_all_regions_adj$dist <- factor(all_taus_dists_emp_cov_all_regions_adj$dist, 
                                                  labels = c(gnorm = TeX('Gen.\\ Normal'),
                                                             norm = 'Normal', t = TeX('$t_3$')))


all_taus_dists_emp_cov_all_regions_adj$var_case <- factor(all_taus_dists_emp_cov_all_regions_adj$var_case, 
                                                      levels = c('Homoscedastic', 'Heteroscedastic'))

ggplot(all_taus_dists_emp_cov_all_regions_adj) + 
  geom_point(aes(x = nom_cov, y = emp_cov)) + 
  facet_grid(var_case ~ dist, labeller = label_parsed) +
  geom_abline(intercept = 0, slope = 1) + 
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(text = element_text(size = 20)) +
  labs(x = 'Nominal Coverage', y = 'Empirical Coverage')
ggsave(paste0('output/sim_data/all_regions_coverage_adj.pdf'),
       height = 8, width = 16, units = 'in', dpi = 600)



#####################
# Region 1
#####################

##########
# Determined by true x
##########
all_taus_dists_emp_cov_truex_region1$dist <- factor(all_taus_dists_emp_cov_truex_region1$dist, 
                                                  labels = c(gnorm = TeX('Gen.\\ Normal'),
                                                             norm = 'Normal', t = TeX('$t_3$')))


all_taus_dists_emp_cov_truex_region1$var_case <- factor(all_taus_dists_emp_cov_truex_region1$var_case, 
                                                      levels = c('Homoscedastic', 'Heteroscedastic'))

ggplot(all_taus_dists_emp_cov_truex_region1) +
  geom_point(aes(x = nom_cov, y = emp_cov)) + 
  facet_grid(var_case ~ dist, labeller = label_parsed) +
  geom_abline(intercept = 0, slope = 1) + 
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(text = element_text(size = 20)) +
  labs(x = 'Nominal Coverage', y = 'Empirical Coverage')
ggsave(paste0('output/sim_data/truex_region1_coverage.pdf'),
       height = 8, width = 16, units = 'in', dpi = 600)



##########
# Determined by estimated x
##########
all_taus_dists_emp_cov_estx_region1$dist <- factor(all_taus_dists_emp_cov_estx_region1$dist, 
                                                    labels = c(gnorm = TeX('Gen.\\ Normal'),
                                                               norm = 'Normal', t = TeX('$t_3$')))


all_taus_dists_emp_cov_estx_region1$var_case <- factor(all_taus_dists_emp_cov_estx_region1$var_case, 
                                                        levels = c('Homoscedastic', 'Heteroscedastic'))

ggplot(all_taus_dists_emp_cov_estx_region1) +
  geom_point(aes(x = nom_cov, y = emp_cov)) + 
  facet_grid(var_case ~ dist, labeller = label_parsed) +
  geom_abline(intercept = 0, slope = 1) + 
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(text = element_text(size = 20)) +
  labs(x = 'Nominal Coverage', y = 'Empirical Coverage')
ggsave(paste0('output/sim_data/estx_region1_coverage.pdf'),
       height = 8, width = 16, units = 'in', dpi = 600)





##########
# Adjusted p-value, Determined by estimated x
##########
all_taus_dists_emp_cov_estx_region1_adj$dist <- factor(all_taus_dists_emp_cov_estx_region1_adj$dist, 
                                                   labels = c(gnorm = TeX('Gen.\\ Normal'),
                                                              norm = 'Normal', t = TeX('$t_3$')))


all_taus_dists_emp_cov_estx_region1_adj$var_case <- factor(all_taus_dists_emp_cov_estx_region1_adj$var_case, 
                                                       levels = c('Homoscedastic', 'Heteroscedastic'))

ggplot(all_taus_dists_emp_cov_estx_region1_adj) +
  geom_point(aes(x = nom_cov, y = emp_cov)) + 
  facet_grid(var_case ~ dist, labeller = label_parsed) +
  geom_abline(intercept = 0, slope = 1) + 
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(text = element_text(size = 20)) +
  labs(x = 'Nominal Coverage', y = 'Empirical Coverage')
ggsave(paste0('output/sim_data/estx_region1_coverage_adj.pdf'),
       height = 8, width = 16, units = 'in', dpi = 600)




#####################
# Region 2
#####################

##########
# Determined by true x
##########
all_taus_dists_emp_cov_truex_region2$dist <- factor(all_taus_dists_emp_cov_truex_region2$dist, 
                                              labels = c(gnorm = TeX('Gen.\\ Normal'),
                                                         norm = 'Normal', t = TeX('$t_3$')))


all_taus_dists_emp_cov_truex_region2$var_case <- factor(all_taus_dists_emp_cov_truex_region2$var_case, 
                                                  levels = c('Homoscedastic', 'Heteroscedastic'))

ggplot(all_taus_dists_emp_cov_truex_region2) +
  geom_point(aes(x = nom_cov, y = emp_cov)) + 
  facet_grid(var_case ~ dist, labeller = label_parsed) +
  geom_abline(intercept = 0, slope = 1) + 
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(text = element_text(size = 20)) +
  labs(x = 'Nominal Coverage', y = 'Empirical Coverage')
ggsave(paste0('output/sim_data/truex_region2_coverage.pdf'),
       height = 8, width = 16, units = 'in', dpi = 600)



##########
# Determined by estimated x
##########
all_taus_dists_emp_cov_estx_region2$dist <- factor(all_taus_dists_emp_cov_estx_region2$dist, 
                                                    labels = c(gnorm = TeX('Gen.\\ Normal'),
                                                               norm = 'Normal', t = TeX('$t_3$')))


all_taus_dists_emp_cov_estx_region2$var_case <- factor(all_taus_dists_emp_cov_estx_region2$var_case, 
                                                        levels = c('Homoscedastic', 'Heteroscedastic'))

ggplot(all_taus_dists_emp_cov_estx_region2) +
  geom_point(aes(x = nom_cov, y = emp_cov)) + 
  facet_grid(var_case ~ dist, labeller = label_parsed) +
  geom_abline(intercept = 0, slope = 1) + 
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(text = element_text(size = 20)) +
  labs(x = 'Nominal Coverage', y = 'Empirical Coverage')
ggsave(paste0('output/sim_data/estx_region2_coverage.pdf'),
       height = 8, width = 16, units = 'in', dpi = 600)



##########
# Adjusted p-value, determined by estimated x
##########
all_taus_dists_emp_cov_estx_region2_adj$dist <- factor(all_taus_dists_emp_cov_estx_region2_adj$dist, 
                                                   labels = c(gnorm = TeX('Gen.\\ Normal'),
                                                              norm = 'Normal', t = TeX('$t_3$')))


all_taus_dists_emp_cov_estx_region2_adj$var_case <- factor(all_taus_dists_emp_cov_estx_region2_adj$var_case, 
                                                       levels = c('Homoscedastic', 'Heteroscedastic'))

ggplot(all_taus_dists_emp_cov_estx_region2_adj) +
  geom_point(aes(x = nom_cov, y = emp_cov)) + 
  facet_grid(var_case ~ dist, labeller = label_parsed) +
  geom_abline(intercept = 0, slope = 1) + 
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(text = element_text(size = 20)) +
  labs(x = 'Nominal Coverage', y = 'Empirical Coverage')
ggsave(paste0('output/sim_data/estx_region2_coverage_adj.pdf'),
       height = 8, width = 16, units = 'in', dpi = 600)




#####################
# Region 3
#####################

##########
# Determined by true x
##########
all_taus_dists_emp_cov_truex_region3$dist <- factor(all_taus_dists_emp_cov_truex_region3$dist, 
                                              labels = c(gnorm = TeX('Gen.\\ Normal'),
                                                         norm = 'Normal', t = TeX('$t_3$')))


all_taus_dists_emp_cov_truex_region3$var_case <- factor(all_taus_dists_emp_cov_truex_region3$var_case, 
                                                  levels = c('Homoscedastic', 'Heteroscedastic'))

ggplot(all_taus_dists_emp_cov_truex_region3) +
  geom_point(aes(x = nom_cov, y = emp_cov)) + 
  facet_grid(var_case ~ dist, labeller = label_parsed) +
  geom_abline(intercept = 0, slope = 1) + 
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(text = element_text(size = 20)) +
  labs(x = 'Nominal Coverage', y = 'Empirical Coverage')
ggsave(paste0('output/sim_data/truex_region3_coverage.pdf'),
       height = 8, width = 16, units = 'in', dpi = 600)


##########
# Determined by estimated x
##########
all_taus_dists_emp_cov_estx_region3$dist <- factor(all_taus_dists_emp_cov_estx_region3$dist, 
                                                    labels = c(gnorm = TeX('Gen.\\ Normal'),
                                                               norm = 'Normal', t = TeX('$t_3$')))


all_taus_dists_emp_cov_estx_region3$var_case <- factor(all_taus_dists_emp_cov_estx_region3$var_case, 
                                                        levels = c('Homoscedastic', 'Heteroscedastic'))

ggplot(all_taus_dists_emp_cov_estx_region3) +
  geom_point(aes(x = nom_cov, y = emp_cov)) + 
  facet_grid(var_case ~ dist, labeller = label_parsed) +
  geom_abline(intercept = 0, slope = 1) + 
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(text = element_text(size = 20)) +
  labs(x = 'Nominal Coverage', y = 'Empirical Coverage')
ggsave(paste0('output/sim_data/estx_region3_coverage.pdf'),
       height = 8, width = 16, units = 'in', dpi = 600)








##########
# Adjusted p-value, determined by estimated x
##########
all_taus_dists_emp_cov_estx_region3_adj$dist <- factor(all_taus_dists_emp_cov_estx_region3_adj$dist, 
                                                   labels = c(gnorm = TeX('Gen.\\ Normal'),
                                                              norm = 'Normal', t = TeX('$t_3$')))


all_taus_dists_emp_cov_estx_region3_adj$var_case <- factor(all_taus_dists_emp_cov_estx_region3_adj$var_case, 
                                                       levels = c('Homoscedastic', 'Heteroscedastic'))

ggplot(all_taus_dists_emp_cov_estx_region3_adj) +
  geom_point(aes(x = nom_cov, y = emp_cov)) + 
  facet_grid(var_case ~ dist, labeller = label_parsed) +
  geom_abline(intercept = 0, slope = 1) + 
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(text = element_text(size = 20)) +
  labs(x = 'Nominal Coverage', y = 'Empirical Coverage')
ggsave(paste0('output/sim_data/estx_region3_coverage_adj.pdf'),
       height = 8, width = 16, units = 'in', dpi = 600)




