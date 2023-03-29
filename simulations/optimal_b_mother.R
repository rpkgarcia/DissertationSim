

# -------------------------------------------------------------------------
# Load Functions ----------------------------------------------------------
# -------------------------------------------------------------------------
url <- "https://raw.githubusercontent.com/rpkgarcia/DissertationSim/main/fixed_b_asymptotics/Fitted_fixed_b/"
source(paste(url, "fitted_cv.R", sep = ""))

url <- "https://raw.githubusercontent.com/rpkgarcia/DissertationSim/main/support_functions/"
source(paste(url, "kernel_stats.R", sep = ""))

# -------------------------------------------------------------------------
# Type 2 Error ------------------------------------------------------------
# -------------------------------------------------------------------------



type2_error <- function(b, rho=0.7, alpha=0.05, big_T=200, d=1, ncp = 5.67, 
                        the_kernel="bartlett", mother = 2){
  the_cv <- qchisq(1-alpha, df = d)
  c2 <- all_c2[[the_kernel]][[mother]] 
  q <- all_q[[the_kernel]]
  c_b <- 2*rho*rho^floor(b*big_T)/(1+rho)
  w_q <- all_w_q(q, rho)
  g_q <- all_g_q[[the_kernel]]
  c<- 1
  
  if(mother == 2){
    c <- 0 
  }
  
  
  p1 <- pchisq(the_cv, df = d, ncp = ncp)
  #p1 <- 0
  p2 <- b*dchisq(the_cv, df = d, ncp = ncp)*the_cv*c2*(the_cv+2)/2      # b to be bigger 
  p3 <- dchisq(the_cv, df = d, ncp = ncp)*the_cv*c_b                    # b to be smaller 
  p4 <-  c*dchisq(the_cv, df = d, ncp = ncp)*the_cv*g_q*w_q*(b*big_T)^(-q) # b to be smaller 
  
  p <- p1 + p2 + p3 + p4
  
  return(p)
}

type2_error(.08)

try_b <- seq(0.01, 1, length.out = 200)
the_big_T <- seq(100, 1000, by =100)
plot(the_big_T, type2_error(b =0.09, big_T = the_big_T), ylab = c("Type 2 Error"))
abline(h = c(0,1), col = "red")


try_b <- seq(0.01, 1, length.out = 200)
plot(try_b, type2_error(b =try_b), ylab = c("Type 2 Error"))
abline(h = c(0,1), col = "red")
# found by optim




# found by me 
# plot power curve 
try_ncp <- seq(.1, 10, length.out = 200)
plot(try_ncp, 1-type2_error(b =0.09, ncp = try_ncp), ylab = c("Type 2 Error"))
abline(h = c(0,1), col = "red")

optim(.015, type2_error, method = "Brent", lower = min(try_b), upper = max(try_b),ncp = 10, 
      mother = 2)


# DOES THE OPTIMAL B depend on the noncentrality parameter???? NO
try_ncp <- 1:6
for(the_ncp in try_ncp){
  try_b <- seq(0.01, 1, length.out = 200)
  
  temp_b <- optim(.015, type2_error, method = "Brent", lower = min(try_b), upper = max(try_b),
                  ncp = the_ncp, mother = 2 )$par
  
  plot(try_b, type2_error(try_b, ncp = the_ncp, mother = 2), ylab = c("Type 2 Error"), 
       main = paste("b=", round(temp_b, 3), "; ncp=", the_ncp), type = "l")
  Sys.sleep(1)
  points(temp_b, type2_error(temp_b, ncp = the_ncp, mother = 2), col = "red")
  abline(v = temp_b, col ="red")
  Sys.sleep(2)
}

# Suns type 2 error -------------------------------------------------------


sun_type2_error <- function(b, rho=0.7, alpha=0.05, big_T=200, d=1, ncp = 5.67, 
                        the_kernel="bartlett", mother = 1, c = 1){
  the_cv <- qchisq(1-alpha, df = d)
  c2 <- all_c2[[the_kernel]][[mother]]
  q <- all_q[[the_kernel]]
  c_b <- 2*rho*rho^floor(b*big_T)/(1+rho)
  w_q <- all_w_q(q, rho)
  g_q <- all_g_q[[the_kernel]]
  
  
  p1 <- pchisq(the_cv, df = d, ncp = ncp)
  #p1 <- 0
  p2 <- ncp*b*dchisq(the_cv, df = (d+2), ncp = ncp)*the_cv*c2/2      # b to be bigger 
  p4 <- -dchisq(the_cv, df = d, ncp = ncp)*the_cv*g_q*w_q*(b*big_T)^(-q)*c # b to be smaller 
  
  p <- p1 + p2+p4
  
  return(p)
}

sun_type2_error(.08, rho =-.7)

# plot power curve 
try_ncp <- seq(.1, 10, length.out = 200)
plot(try_ncp, 1-sun_type2_error(b =0.9, ncp = try_ncp), ylab = c("Type 2 Error"))
abline(h = c(0,1), col = "red")


try_b <- seq(0.01, 1, length.out = 200)
plot(try_b, sun_type2_error(try_b, rho = -.7))
abline(h = c(0,1), col = "red")


the_big_T <- seq(100, 1000, by =100)
plot(the_big_T, sun_type2_error(b =0.09, big_T = the_big_T), ylab = c("Type 2 Error"))
abline(h = c(0,1), col = "red")

optim(.015, sun_type2_error, method = "Brent", lower = min(try_b), upper = max(try_b))



# DOES THE OPTIMAL B depend on the noncentrality parameter???? YES 
try_ncp <- seq(1, 20, length.out = 10)
for(the_ncp in try_ncp){
  try_b <- seq(0.01, 1, length.out = 200)
  
  temp_b <- optim(.015, sun_type2_error, method = "Brent", lower = min(try_b), 
                  upper = max(try_b),ncp = the_ncp)$par
  
  plot(try_b, sun_type2_error(try_b, ncp = the_ncp), ylab = c("Type 2 Error"), 
       main = paste("SUN - b=", round(temp_b, 3), "; ncp=", the_ncp), type = "l")
  Sys.sleep(1)
  points(temp_b, type2_error(temp_b, ncp = the_ncp), col = "red")
  abline(v = 0.085, col ="red")
  Sys.sleep(2)
}

# both --------------------------------------------------------------------
try_b <- seq(0.01, 1, length.out = 200)
plot(try_b, sun_type2_error(try_b), type = "l", ylim = c(0,1))                         # Suns - Mother 
lines(try_b, sun_type2_error(try_b, c = 0), type = "l", ylim = c(0,1), col = "orange") # Suns - Zero Lugsail
lines(try_b, type2_error(try_b), col = "blue")              # AR(1) - Mother 
lines(try_b, type2_error(try_b, mother = 2), col = "green") # AR(1) - Zero Lugsail 
abline(h = c(0,1), col = "red")


b_opt <- optim(.015, type2_error, method = "Brent", lower = min(try_b), upper = max(try_b))$par
b_opt_lug<- optim(.015, type2_error, method = "Brent", lower = min(try_b), upper = max(try_b), mother = 2)$par
sun_b_opt <- optim(.015, sun_type2_error, method = "Brent", lower = min(try_b), upper = max(try_b))$par

points(b_opt, type2_error(b_opt))
points(b_opt_lug, type2_error(b_opt_lug, mother = 2))
points(sun_b_opt, sun_type2_error(sun_b_opt))

c(b_opt, sun_b_opt, b_opt_lug)
