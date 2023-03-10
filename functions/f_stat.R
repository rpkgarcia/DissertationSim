# -------------------------------------------------------------------------
# Calculate F-Statistics --------------------------------------------------
# -------------------------------------------------------------------------

# Calculate the F test statistic. 

# b: the bandwidth, proportion of estimated autocovs that have non-zero weight
# the_sim_data: 
#     - Contains `big_T` simulated dependent random vectors of dimension (d x 1).
#     - Should already be centered with hypothesized or estimated means.
# the_means: 
#     - The estimated mean vector. 
#     - Should be of length d.
#     - Need this because the_sim_data is already centered. 
# all_autocovariances: 
#     - Rows are correspond to estimated autovoc at lag [0, ..., (big_T -1)]
#     - Columns correspond to the vectorization of the estimated autocov matrix: 
#          R11, R12, R13, ..., R1d, R21, R22, ..., R2d, ..., Rd1, ...Rdd
# kernel: the name of the kernel function 
# lugsail_parameters: 
#     - a named list, list(r = 1, c= 0) that contains the lugsail parameters
#     - default is non-lugsail  

F_stats <- function(b, sim_data, the_means, all_autocovariances, 
                   kernel, lugsail_parameters = list(r = 1, c= 0)){
  
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
    F_stat = the_means[1:d]%*% omega_hat_inv %*%the_means[1:d]/d
  })
  
  return(F_stat)
}
