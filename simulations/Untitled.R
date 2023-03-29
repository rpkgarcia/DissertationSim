# Compare optimal bs


rho=0.7
alpha=0.05
big_T=200
kernel=1  # Barlett
mother=1 
d=1


optimal_b_zero(rho=rho, alpha=alpha, big_T=big_T, kernel=1, mother=2, d=d)
optim(.015, type2_error, method = "Brent", lower = .001, upper = 1, 
      rho=rho, alpha=alpha, big_T=big_T, d=d, ncp = 5.67, the_kernel="bartlett", mother = 1)$par
optim(.015, type2_error, method = "Brent", lower = .001, upper = 1, 
      rho=rho, alpha=alpha, big_T=big_T, d=d, ncp = 5.67, the_kernel="bartlett", mother = 2)$par




try_b <- seq(0.01, 1, length.out = 200)
plot(try_b, type2_error(b =try_b, mother = 2), ylab = c("Type 2 Error"))
abline(h = c(0,1), col = "red")
# found by optim
opt_b_optim <- optim(.1, type2_error, method = "Brent", lower = .001, upper = 1, 
                     rho=rho, alpha=alpha, big_T=big_T, d=d, ncp = 5.67, the_kernel="bartlett", mother = 2)$par
abline(v = opt_b_optim, col = "red")

# found by me 
opt_b <- optimal_b_zero(rho=rho, alpha=alpha, big_T=big_T, kernel=1, mother=2, d=d)
abline(v = opt_b, col = "blue")



# plot the bandwidth as a function of T ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
the_big_T <- c(seq(100, 900, by = 100 ), seq(1000, 10000, by =1000))
opt_b_optim <- sapply(the_big_T, function(the_T){
                    optim(.015, type2_error, method = "Brent", lower = .001, upper = 1, 
                    rho=rho, alpha=alpha, big_T=the_T, d=d, ncp = 5.67, the_kernel="bartlett", mother = 2)$par})

opt_b_mother_optim <- sapply(the_big_T, function(the_T){
  optim(.015, sun_type2_error, method = "Brent", lower = .001, upper = 1, 
        rho=rho, alpha=alpha, big_T=the_T, d=d, ncp = 5.67, the_kernel="bartlett", mother = 1)$par})

opt_b_zero <- sapply(the_big_T, function(the_T){
  optimal_b_zero(rho=rho, alpha=alpha, big_T=the_T, kernel=1, mother=2, d=d)
})

plot(the_big_T, opt_b_optim, col = "black", main = "Optimal b, function of T", 
     xlab = "T", ylab = "b", ylim = range(opt_b_optim, opt_b_mother_optim, opt_b_zero), 
     lty = 2, type = "l")
lines(the_big_T, opt_b_zero, col = "blue")
lines(the_big_T, opt_b_mother_optim, col = "red")

legend("topright", 
       c("Mother-Optim",
         "Zero-Optim", 
         "Zero-Analytical"), 
       col = c("red", "black", "blue"), 
       lty = c(1, 2, 1))

