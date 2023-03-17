# Compare optimal bs


rho=0.7
alpha=0.05
big_T=200
kernel=1  # Barlett
mother=1 
d=1

find_detla_sq <- function(ncp, rho=0.7, alpha=0.05, big_T=200, d=1){
  b_opt <- optimal_b_zero(rho=rho, alpha=alpha, big_T=big_T, kernel=1, mother=2, d=d)
  optim_b_opt <- optim(0.01, type2_error, method = "Brent", lower = .001, upper = 1, 
                       rho=rho, alpha=alpha, big_T=big_T, d=d, ncp = ncp, 
                       the_kernel="bartlett", mother = 2)$par
  return(abs(b_opt - optim_b_opt))
}

optim(5.67, find_delta_sq, method ="Brent", lower = 0.001, upper = 10)

optimal_b_zero(rho=rho, alpha=alpha, big_T=big_T, kernel=1, mother=2, d=d)
optim(.015, type2_error, method = "Brent", lower = .001, upper = 1, 
      rho=rho, alpha=alpha, big_T=big_T, d=d, ncp = 5.67, the_kernel="bartlett", mother = 1)$par
optim(.015, type2_error, method = "Brent", lower = .001, upper = 1, 
      rho=rho, alpha=alpha, big_T=big_T, d=d, ncp = 5.67, the_kernel="bartlett", mother = 2)$par
