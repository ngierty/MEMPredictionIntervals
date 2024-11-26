##########################################################################################
# Create output for paper from mass-radius relationship analysis in MR_02_calc_ncss.R
##########################################################################################

##########################################
# Libraries, functions, and global variables
##########################################

library(ggplot2)
library(reshape2)


# Set the break points based on Ma (2019) and Ma (2021)
zetas <- c(1.5, 3.84, 8, 13)
log10_zetas <- log10(zetas)


##########################################
# Load the non-conformity scores and calculate regional-based p-values
##########################################

test_idx <- 1
ntest <- 300
ncs_region_pvalue <- rep(NA, ntest)
test_xhats <- rep(NA, ntest)
test_xhat_region <- rep(NA, ntest)

for(test_idx in 1:ntest){
  # Load in the non-conformity scores and estimated x-values
  load(file = paste0('inter/mr_data/calc_ncss_vals__test_idx=', test_idx, '.RData'))
  
  # Extract the non-conformity scores
  ncss <- unlist(calc_ncss_vals[1,])
  n_p_1 <- length(ncss)
  
  # Extract the estimated x-values
  xhats <- unlist(calc_ncss_vals[2,])
  
  # Save the test x values
  test_xhats[test_idx] <- xhats[n_p_1]
  
  
  # Calculate the p-value of the non-conformity score based on the estimated region
  estx_region1_flag <- (xhats < log10_zetas[1])
  estx_region2_flag <- (log10_zetas[1] <= xhats) & (xhats < log10_zetas[2])
  estx_region3_flag <- (log10_zetas[2] <= xhats) & (xhats < log10_zetas[3])
  estx_region4_flag <- (log10_zetas[3] <= xhats) & (xhats < log10_zetas[4])
  estx_region5_flag <- (xhats >= log10_zetas[4])
  
  if(estx_region1_flag[n_p_1]){
    
    test_xhat_region[test_idx] <- 1
    ncs_region_pvalue[test_idx] <- mean(ncss[estx_region1_flag] > ncss[n_p_1])
    
  }else if(estx_region2_flag[n_p_1]){
    
    test_xhat_region[test_idx] <- 2
    ncs_region_pvalue[test_idx] <- mean(ncss[estx_region2_flag] > ncss[n_p_1])
    
  }else if(estx_region3_flag[n_p_1]){
    
    test_xhat_region[test_idx] <- 3
    ncs_region_pvalue[test_idx] <- mean(ncss[estx_region3_flag] > ncss[n_p_1])
    
  }else if(estx_region4_flag[n_p_1]){
    
    test_xhat_region[test_idx] <- 4
    ncs_region_pvalue[test_idx] <- mean(ncss[estx_region4_flag] > ncss[n_p_1])
    
  }else{
    
    test_xhat_region[test_idx] <- 5
    ncs_region_pvalue[test_idx] <- mean(ncss[estx_region5_flag] > ncss[n_p_1])
  }
}

# Save into a dataframe
pvalues_xhats <- data.frame(test_point = 1:ntest, xhat = test_xhats,
                            region = test_xhat_region,
                            pval = ncs_region_pvalue)

# Calculate the 95% empirical coverage for each region
  # not enough data points to split into different nominal coverages as below
region_95_emp_cov <- rep(NA, 5)
for(i in 1:5){
  
  temp <- pvalues_xhats[pvalues_xhats$region == i,]
  region_95_emp_cov[i] <- mean(temp$pval >= 0.05)
}

region_95_emp_cov <- data.frame(region = 1:5, emp_cov = region_95_emp_cov)
write.csv(region_95_emp_cov, 'output/mr_data/region_95_emp_cov.csv',
          row.names = FALSE)

##########################################
# Plot the nominal vs empirical coverages for all regions
  # for a small set of nominal coverages
##########################################

# Choose a smaller region of type 1 errors/nominal coverages to consider
type1_err <- seq(0.01, 0.2, length.out = 10)
nom_cov <- 1 - type1_err
nerr <- length(type1_err)

# Randomly split the test points (and corresponding p-values) into ones 
  # that will be used for each coverage level
set.seed(7958)
err_split <- data.frame(fold = rep(1:nerr, each = ntest/nerr))
err_split$test_point <- NA
unselected_test_points <- 1:ntest
for(i in 1:nerr){
  sample_points <- sample(unselected_test_points, ntest/nerr, replace = FALSE)
  err_split[err_split$fold == i, 'test_point'] <- sample_points
  unselected_test_points <- setdiff(unselected_test_points, sample_points)
}

# Split p-values into different folds for calculating nominal coverages
pvalues_split <- merge(pvalues_xhats, err_split)
pvalues_split <- pvalues_split[order(pvalues_split$fold),]
pvalues_split$timevar <- rep(1:(ntest/nerr), nerr)


  # drops the test index
pvalues_split <- pvalues_split[,c('fold', 'timevar', 'pval')]
pvalues_split_wide <- reshape(pvalues_split, direction = 'wide',
                              idvar = 'fold', timevar = 'timevar')

pvalues_split_wide <- subset(pvalues_split_wide, select=-c(fold))
all_regions_cov <- data.frame(nom_cov)
all_regions_cov$emp_cov <- NA

for(i in 1:nerr){
  # All regions
  all_regions_cov$emp_cov[i] <- mean(pvalues_split_wide[i,] > type1_err[i], na.rm = TRUE)
}

ymin <- min(all_regions_cov$emp_cov, 0.8)

ggplot(all_regions_cov) + 
  geom_point(aes(x = nom_cov, y = emp_cov)) + 
  geom_abline(intercept = 0, slope = 1) + 
  theme(text = element_text(size = 30)) +
  ylim(c(ymin, 1)) +
  labs(x = 'Nominal Coverage', y = 'Empirical Coverage')
ggsave(paste0('output/mr_data/all_regions_coverage.pdf'),
       height = 8, width = 10, units = 'in', dpi = 600)

















