# -------------------------------------------------------------------------
# Calculate F-Statistics --------------------------------------------------
# -------------------------------------------------------------------------

# Calculate the F test statistic for dimensions 1, ..., d_max
# Under the 

# b: the bandwidth, proportion of estimated autocovs that have non-zero weight
# the_sim_data: 
#     - Contains `big_T` simulated dependent random vectors of dimension (d_max).
#     - Should already be centered with hypothesized or estimated means.
# the_means: 
#     - The estimated mean vector.  
#     - Should be of length d.
#     - Need this because the_sim_data is already centered. 
# all_autocovariances: 
#     - Rows are correspond to estimated autocov (R) at lag [0, ..., (big_T -1)]
#     - Columns correspond to the vectorization of the estimated autocov matrix: 
#          R11, R12, R13, ..., R1d, R21, R22, ..., R2d, ..., Rd1, ...Rdd
# kernel: the name of the kernel function 
# lugsail_parameters: 
#     - a named list, list(r = 1, c= 0) that contains the lugsail parameters
#     - default is non-lugsail  

F_stats <- function(b, sim_data, the_means, all_autocovariances, 
                    kernel, lugsail_parameters = list(r = 1, c= 0),
                    null_means = rep(0, ncol(sim_data))){
  
  d_max <- ncol(sim_data)
  big_T <- nrow(sim_data)
  M <- b*nrow(sim_data)
  
  # Wieghts matrix 
  W <- rep(0, nrow(all_autocovariances))
  nonzero_weights <- sapply(0:(M)/(M+1), lugsail, kernel = kernel, 
                            lugsail_parameters = lugsail_parameters)
  
  # Everything beyond M has a zero weight
  W[1:length(nonzero_weights)] <- nonzero_weights
  
  # Calculate omega_hat 
  omega_hat <- matrix(colSums(W * all_autocovariances)/big_T, d_max, d_max) 
  
  # Calculate F-statistics for d = 1, ...., d_max for each simulation 
  F_stat <- sapply(1:d_max, function(d){
    omega_hat_inv = solve(omega_hat[1:d, 1:d])
    num <- (the_means- null_means)[1:d]
    F_stat = (num%*% omega_hat_inv %*%num)/d
  })
  
  return(F_stat)
}



# -------------------------------------------------------------------------
# Main --------------------------------------------------------------------
# -------------------------------------------------------------------------


# Generates a null data set.  Calculates the autocovariance matrices. 
# Calculates the F-statistic for dimensions (1,..., d_max). 
# Need to repeat this simulation many times to then find what asymptotic CV is. 


simulate_f_stat <- function(big_T = 1000, d_max = 12){
  
  # Simulate the data 
  sim_data <- matrix(rnorm(d_max*big_T), nrow = big_T, ncol = d_max)
  the_means <- colMeans(sim_data)
  sim_data <- apply(sim_data, 1, function(row) row - the_means)
  sim_data <- t(sim_data)
  
  
  # ------- AutoCovariance Matrices  -------
  # Each row [#, ] is the lag (0, ..., big_T-1)
  # Each column [ , #] is a component of the vectorized autocov matrix. 
  #    R11, R12,  ..., R1d; R21, R22, ..., R2d; ...; Rd1, ...Rdd
  all_autocovariances <-sapply(0:(big_T-1), R, the_sim_data = sim_data)
  all_autocovariances <- t(all_autocovariances)

  
  
  # ------- F-statistics for various b values (an 2-dimensional array) -------
  # Each row    [ #,  ] is a different b value
  # Each column [  , #] is a different m
  # colnames(F_stats) = paste("m", 1:ncol(sim_data), sep = "")
  # rownames(F_stats) = paste("b=", try_b, sep="")
  
  # ------- BARTLETT ------- 
  F_stats <- sapply(try_b, F_stats, sim_data=sim_data, 
                   the_mean = the_means, 
                   kernel = bartlett, 
                   all_autocovariances = all_autocovariances, 
                   lugsail_parameters = list(r = 1, c = 0))
  F_stats_bartlett <- t(F_stats)
  
  F_stats <- sapply(try_b, F_stats, sim_data=sim_data, 
                   the_mean = the_means, 
                   kernel = bartlett, 
                   all_autocovariances = all_autocovariances, 
                   lugsail_parameters = list(r = 2, c = .5))
  F_stats_bartlett_lug <- t(F_stats)
  
  # #  ------- PARZEN ------- 
  # F_stats = sapply(try_b, F_stats, sim_data=sim_data, 
  #                  the_mean = the_means, 
  #                  kernel = parzen, 
  #                  all_autocovariances = all_autocovariances, 
  #                  lugsail_parameters = F)
  # F_stats_parzen = t(F_stats)
  # 
  # F_stats = sapply(try_b, F_stats, sim_data=sim_data, 
  #                  the_mean = the_means, 
  #                  kernel = parzen, 
  #                  all_autocovariances = all_autocovariances, 
  #                  lugsail_parameters = lugsail_parameters)
  # F_stats_parzen_lug = t(F_stats)
  # 
  # # -------  QS ------- 
  # F_stats = sapply(try_b, F_stats, sim_data=sim_data, 
  #                  the_mean = the_means, 
  #                  kernel = qs, 
  #                  all_autocovariances = all_autocovariances, 
  #                  lugsail_parameters = F)
  # F_stats_qs = t(F_stats)
  # 
  # F_stats = sapply(try_b, F_stats, sim_data=sim_data, 
  #                  the_mean = the_means, 
  #                  kernel = qs, 
  #                  all_autocovariances = all_autocovariances, 
  #                  lugsail_parameters = lugsail_parameters)
  # F_stats_qs_lug = t(F_stats)
  
  
  
  # ------- Store everything -------
  # [ #,  ,  ] is a different b value
  # [  , #,  ] is a different d
  # [  ,  , #] is a kernel 
  F_stats_all <- array(c(F_stats_bartlett, F_stats_bartlett_lug), 
                      # F_stats_parzen, F_stats_parzen_lug, 
                      # F_stats_qs, F_stats_qs_lug), 
                      dim=c(length(try_b), d_max, 2))
  
  dimnames(F_stats_all)[[3]] <- c("Bartlett", "Bartlett Lugsail")  
  # "Parzen", "Parzen Lugsail", 
  # "QS", "QS Lugsail")
  return(F_stats_all)
}

# -------------------------------------------------------------------------
# Simulation Parameters (Adjust me) ---------------------------------------
# -------------------------------------------------------------------------

#  What bs to use 
try_b <-  seq(0.005, .99, by = 0.005)

# How many replicates
# KV005 used 50,0000
num_replicates <- 100

# Sample size of each replicate
big_T = 1000

# Maximum number of dimensions
mdmax = 12


# -------------------------------------------------------------------------
# Run Simulation  ---------------------------------------------------------
# -------------------------------------------------------------------------
all_F_stats = replicate(num_replicates, 
                        simulate_f_stat(big_T = big_T, m_max = m_max))

# [ #,  ,  ,  ] : rows, the different b values [1:length(try_b)]
# [  , #,  ,  ] : columns, the different d values [1:max_d]
# [  ,  , #,  ] : the original kernel (1), and lugsail (2) [1:2]
# [  ,  ,  , #] : each simulation [1:num_replicates]

dimnames(all_F_stats)[[1]] = paste("b=", try_b, sep="")
dimnames(all_F_stats)[[2]] = paste("m", 1:m_max, sep="")
dimnames(all_F_stats)[[3]] = c("Bartlett", "Bartlett_Lugsail") 
# "Parzen", "Parzen_Lugsail", 
# "QS", "QS_Lugsail")
dimnames(all_F_stats)[[4]] = paste("sim", 1:num_replicates, sep="")

# -------------------------------------------------------------------------
# Save Results ------------------------------------------------------------
# -------------------------------------------------------------------------

save_results <- function(kernel){
  setwd("Fixed_b_CV_Tables")
  bartlett_F_stats <- all_F_stats[ , , kernel, ]
  
  t90 <- apply(bartlett_F_stats, 1:2, quantile, probs = 0.90)
  write.csv(t90, paste(kernel, "_10.csv", sep = ""))
  
  t95 <- apply(bartlett_F_stats, 1:2, quantile, probs = 0.95)
  write.csv(t95, paste(kernel, "_05.csv", sep = ""))
  
  t97.5 <- apply(bartlett_F_stats, 1:2, quantile, probs = 0.975)
  write.csv(t97.5, paste(kernel, "_025.csv", sep = ""))
  
  t99 <- apply(bartlett_F_stats, 1:2, quantile, probs = 0.99)
  write.csv(t99, paste(kernel, "_01.csv", sep = ""))
  setwd("..")
}

#  this is how I get all the critical values 
sapply(dimnames(all_F_stats)[[3]], save_results)

