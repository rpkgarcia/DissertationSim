# plot_density 

plot_b <- seq(.001, .25, length.out = 200)

alpha <- .05

d <- 1


density_mother <- get_cv(d, alpha, "Bartlett", FALSE, plot_b)
density_zero   <- get_cv(d, alpha, "Bartlett", TRUE, plot_b)
density_smallb <- get_cv(d, alpha, "Bartlett", FALSE, 0)

plot(plot_b, density_mother, col = "red", xlab = "b", ylab = "density", type = "l", lty = 2)
lines(plot_b, density_zero, col = "blue")



Bartlett_05 <- read.csv("Bartlett_d1.csv")
Barltett_Lugsail_05 <- read.csv("Bartlett_Lugsail_d1.csv")

density_Bartlett_05 <- density(Bartlett_05[,19], from = 0 , to = 100)  #optimal b for mother 
density_Bartlett_Lugsail_05 <- density(Barltett_Lugsail_05[,9], from=0, to = 100) #optimal b for zero
plot(density_Bartlett_05, col = "blue", ylab = "density", xlab = "X", 
     main = "density", xlim = c(0, 6), lty = 2, lwd = 3)
lines(density_Bartlett_Lugsail_05, col = "red", lty = 3, lwd = 3)
curve(dchisq(x, 1), add = TRUE)

# bart_cv <- get_cv(d, alpha, "Bartlett", F, 0.09)
# bart_lug_cv <- get_cv(d, alpha, "Bartlett", T, 0.04)
bart_cv <- 5.426733
bart_lug_cv <- 4.39494


abline(v = bart_cv, col = "blue")
abline(v = bart_lug_cv, col = "red")
abline(v = qchisq(1-alpha, d))

legend("topright", 
       c("Chi-sq", 
         "Bartlett-Mother", 
         "Bartlett-Zero"), 
       col = c("Black", "blue", "red"), 
       lty = c(1, 2, 3), 
       lwd = c(1, 3, 3))

