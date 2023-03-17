# -------------------------------------------------------------------------
# Load Functions ----------------------------------------------------------
# -------------------------------------------------------------------------
url <- "https://raw.githubusercontent.com/rpkgarcia/DissertationSim/main/support_functions/"

source(paste(url, "est_autocov.R", sep = ""))
source(paste(url, "kernels.R", sep = ""))

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
                    the_kernel, lugsail_parameters = list(r = 1, c= 0),
                    null_means = rep(0, ncol(sim_data))){
  
  d_max <- ncol(sim_data)
  big_T <- nrow(sim_data)
  M <- b*nrow(sim_data)
  
  # Wieghts matrix 
  W <- rep(0, nrow(all_autocovariances))
  nonzero_weights <- sapply(0:(M)/(M+1), lugsail, the_kernel = the_kernel, 
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
  the_means <- colMeans(the_sim_data)
  the_sim_data <- apply(the_sim_data, 1, function(row) row - the_means)
  the_sim_data <- t(the_sim_data)
  
  # ------- AutoCovariance Matrices  -------
  # [#, ] the lag (0, ..., big_T-1)
  # [ , #]  a component of the vectorized autocov matrix. 
  #    R11, R12,  ..., R1d; R21, R22, ..., R2d; ...; Rd1, ...Rdd
  all_autocovariances <-sapply(0:(big_T-1), R, the_sim_data = sim_data)
  all_autocovariances <- t(all_autocovariances)

  
  
  # ------- F-statistics for various b values (an 2-dimensional array) -------
  # [ #,  ] a different b value
  # [  , #] a different d
  
  # ------- BARTLETT ------- 
  F_stats <- sapply(try_b, F_stats, sim_data=sim_data, 
                   the_mean = the_means, 
                   the_kernel = bartlett, 
                   all_autocovariances = all_autocovariances, 
                   lugsail_parameters = list(r = 1, c = 0))
  F_stats_bartlett <- t(F_stats)
  
  F_stats <- sapply(try_b, F_stats, sim_data=sim_data, 
                   the_mean = the_means, 
                   the_kernel = bartlett, 
                   all_autocovariances = all_autocovariances, 
                   lugsail_parameters = list(r = 2, c = .5))
  F_stats_bartlett_lug <- t(F_stats)
  
  # #  ------- PARZEN ------- 
  # F_stats = sapply(try_b, F_stats, sim_data=sim_data, 
  #                  the_mean = the_means, 
  #                  the_kernel = parzen, 
  #                  all_autocovariances = all_autocovariances, 
  #                  lugsail_parameters = F)
  # F_stats_parzen = t(F_stats)
  # 
  # F_stats = sapply(try_b, F_stats, sim_data=sim_data, 
  #                  the_mean = the_means, 
  #                  the_kernel = parzen, 
  #                  all_autocovariances = all_autocovariances, 
  #                  lugsail_parameters = lugsail_parameters)
  # F_stats_parzen_lug = t(F_stats)
  # 
  # # -------  QS ------- 
  # F_stats = sapply(try_b, F_stats, sim_data=sim_data, 
  #                  the_mean = the_means, 
  #                  the_kernel = qs, 
  #                  all_autocovariances = all_autocovariances, 
  #                  lugsail_parameters = F)
  # F_stats_qs = t(F_stats)
  # 
  # F_stats = sapply(try_b, F_stats, sim_data=sim_data, 
  #                  the_mean = the_means, 
  #                  the_kernel = qs, 
  #                  all_autocovariances = all_autocovariances, 
  #                  lugsail_parameters = lugsail_parameters)
  # F_stats_qs_lug = t(F_stats)
  
  
  
  # ------- Store everything -------
  # [ #,  ,  ] a different b value
  # [  , #,  ] a different d
  # [  ,  , #] a the_kernel 
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
num_replicates <- 10000

# Sample size of each replicate
big_T = 1000

# Maximum number of dimensions
# Should do 12
d_max = 2


# -------------------------------------------------------------------------
# Run Simulation  ---------------------------------------------------------
# -------------------------------------------------------------------------
all_F_stats = replicate(num_replicates, 
                        simulate_f_stat(big_T = big_T, d_max = d_max))

# [ #,  ,  ,  ] : rows, the different b values [1:length(try_b)]
# [  , #,  ,  ] : columns, the different d values [1:max_d]
# [  ,  , #,  ] : the original kernel (1), and lugsail (2) [1:2]
# [  ,  ,  , #] : each simulation [1:num_replicates]

dimnames(all_F_stats)[[1]] <- try_b
dimnames(all_F_stats)[[2]] <- 1:d_max
dimnames(all_F_stats)[[3]] <- c("Bartlett", "Bartlett_Lugsail") 
# "Parzen", "Parzen_Lugsail", 
# "QS", "QS_Lugsail")
dimnames(all_F_stats)[[4]] <- paste("sim", 1:num_replicates, sep="")

# -------------------------------------------------------------------------
# Create CV Tables --------------------------------------------------------
# -------------------------------------------------------------------------

save_cv <- function(the_kernel){
  F_stats <- all_F_stats[ , , the_kernel, ]
  
  t90 <- apply(F_stats, 1:2, quantile, probs = 0.90)
  write.csv(t90, paste(the_kernel, "_10.csv", sep = ""))
  
  t95 <- apply(F_stats, 1:2, quantile, probs = 0.95)
  write.csv(t95, paste(the_kernel, "_05.csv", sep = ""))
  
  t97.5 <- apply(F_stats, 1:2, quantile, probs = 0.975)
  write.csv(t97.5, paste(the_kernel, "_025.csv", sep = ""))
  
  t99 <- apply(F_stats, 1:2, quantile, probs = 0.99)
  write.csv(t99, paste(the_kernel, "_01.csv", sep = ""))
}

#  this is how I get all the critical values 
setwd("Fixed_b_CV_Tables")
sapply(dimnames(all_F_stats)[[3]], save_cv)
setwd("..")


# -------------------------------------------------------------------------
# Save for Distribution Plots ---------------------------------------------
# -------------------------------------------------------------------------

# Saves a data frame for a specific the_kernel and dimension. 
# Keeps the F-statistics so I can use them to make distribution plots later. 

save_distr <- function(the_kernel, d){
  F_stats <- t(all_F_stats[ , d, the_kernel, ])
  write.csv(F_stats, paste(the_kernel, "_d", d, ".csv", sep = ""))
}

#  this is how I get all the critical values 
setwd("Fixed_b_distribution")
sapply(1:dim(all_F_stats)[2], function(d){
  sapply(dimnames(all_F_stats)[[3]], save_distr, d = d)
})
setwd("..")
