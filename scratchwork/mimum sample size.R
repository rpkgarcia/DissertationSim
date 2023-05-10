c2 <- 1# Bartlett lugsail 

try_rho <- seq(0.95, 0.999, by= 0.01)
names(try_rho) <- try_rho
try_T <- seq(1, 10, by = 1)
names(try_T) <- try_T
d <- 3

p1 <- delta_sq*dchisq(z, (d+2), ncp = delta_sq)/delta_sq*dchisq(z, d, ncp = delta_sq)


sapply(try_T, function(the_big_T){
  (c2*(1- try_rho)*(delta_sq*p1 - (d-1)/2) + 1)^(1/the_big_T)
}) 


