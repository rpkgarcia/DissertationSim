library(MASS)

VAR1_data <- function(rho_y=.7,  d_max = 12, e_sd=diag(d_max), big_T=200, mean_vec = rep(0, d_max)){
  
  # Generate the data.  Initial value is 0
  sim_data <- matrix(0, nrow = big_T, ncol = d_max)
  
  # The rest of the values
  for(t in 2:big_T){
    e_t <- mvrnorm(1, mu = rep(0, d_max), Sigma = e_sd)
    
    sim_data[t, ] <- mean_vec + rho_y*sim_data[c(t-1), ] + e_t
  }
  
  return(sim_data)
}
