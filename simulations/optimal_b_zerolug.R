
url <- "https://raw.githubusercontent.com/rpkgarcia/DissertationSim/main/support_functions/"
source(paste(url, "kernel_stats.R", sep = ""))

optimal_b_zero <- function(rho=0.7, alpha=0.05, big_T=200, kernel=1, mother=1, d=1){
  the_cv <- qchisq(1-alpha, df = d)
  c2 <- all_c2[[kernel]][[mother]]
  
  a1 <- -c2*(2*the_cv + 1)*(1+ rho)/(8*rho*log(rho)*big_T)
  
  the_b <- log(a1)/(big_T*log(rho))
  return(the_b)
}
