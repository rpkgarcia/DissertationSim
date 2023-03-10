# -------------------------------------------------------------------------
# Estimate Autocovariance Matrix  -----------------------------------------
# -------------------------------------------------------------------------
# Calculates the lag-h autocovariance matrix 

# the_sim_data: 
#     - Contains `big_T` simulated dependent random vectors of dimension (d x 1).
#     - Should already be centered with hypothesized or estimated means.
# h: is the lag-h autocovariance calculated

R <- function(h, the_sim_data){
  big_T <- nrow(the_sim_data)
  index <- 1:(big_T -h)
  
  # Already centered
  est <- lapply(index, function(i){
    est <- the_sim_data[i, ]%*%t(the_sim_data[(i+h), ])
  })
  
  # Sum together
  autocov_s <- Reduce('+', est)/big_T
  
  # Because of symmetry 
  if(h!=0){
    autocov_s <- autocov_s + t(autocov_s)
  }
  
  return(autocov_s)
}
