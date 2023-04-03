url <- "https://raw.githubusercontent.com/rpkgarcia/DissertationSim/main/support_functions/"
source(paste(url, "kernels.R", sep = ""))
source(paste(url, "plot_settings.R", sep = ""))


# -------------------------------------------------------------------------
# Mother  -----------------------------------------------------------------
# -------------------------------------------------------------------------



start_plot("mother_kernels")
par(mfrow = c(1,3))
x <- seq(-1, 1, length.out = 1000)
plot(x, sapply(x, bartlett), 
     ylab = "", 
     xlab = "x", 
     main = "Bartlett", type = "l",
     xaxt="n")
axis(side=1, at=seq(-1, 1, by = 1), labels = T)
axis(side=1, at=seq(-1, 1, by = 0.25), labels = F)

plot(x, sapply(x, parzen), 
     ylab = "", 
     xlab = "x", 
     main = "Parzen", type = "l" ,
     xaxt="n")
axis(side=1, at=seq(-1, 1, by = 1), labels = T)
axis(side=1, at=seq(-1, 1, by = 0.25), labels = F) 
plot(x, sapply(x, qs), 
     ylab = "", 
     xlab = "x", 
     main = "Quadratic Spectral", type = "l",
     xaxt="n")
axis(side=1, at=seq(-1, 1, by = 1), labels = T)
axis(side=1, at=seq(-1, 1, by = 0.25), labels = F)
dev.off()


# -------------------------------------------------------------------------
# All  -----------------------------------------------------------------
# -------------------------------------------------------------------------


start_plot("mother_kernels")
par(mfrow = c(1,3))
x <- seq(0, 1, length.out = 1000)
range_y <- c(0, 1.4)
plot(x, sapply(x, bartlett), 
     ylab = "", 
     xlab = "x", 
     main = "Bartlett", type = "l",
     xaxt="n", 
     ylim = range_y)
lines(x, sapply(x, lugsail, 
                lugsail_parameters = list(r = 2, c = 1/2), 
                the_kernel= bartlett), 
      col = "red", lty = 2)
lugsail_parameters <- get_lugsail_parameters(200, 1, method = "Adaptive")
lines(x, sapply(x, lugsail, 
                lugsail_parameters =lugsail_parameters, 
                the_kernel = bartlett), 
      col = "blue", lty = 3)
lugsail_parameters <- get_lugsail_parameters(big_T = 200, q=1, method = "Over")
lines(x, sapply(x, lugsail, 
                lugsail_parameters =lugsail_parameters, 
                the_kernel = bartlett), 
      col = "green", lty = 4)
axis(side=1, at=seq(0, 1, by = 0.5), labels = T)
axis(side=1, at=seq(0, 1, by = 0.25), labels = F)

plot(x, sapply(x, parzen), 
     ylab = "", 
     xlab = "x", 
     main = "Parzen", type = "l" ,
     xaxt="n",
     ylim = range_y)
lines(x, sapply(x, lugsail, 
                lugsail_parameters = list(r = 2, c = 1/2), 
                the_kernel= parzen), 
      col = "red", lty = 2)
lugsail_parameters <- get_lugsail_parameters(200, 2, method = "Adaptive")
lines(x, sapply(x, lugsail, 
                lugsail_parameters =lugsail_parameters, 
                the_kernel = parzen), 
      col = "blue", lty = 3)
lugsail_parameters <- get_lugsail_parameters(big_T = 200, q=2, method = "Over")
lines(x, sapply(x, lugsail, 
                lugsail_parameters =lugsail_parameters, 
                the_kernel = parzen), 
      col = "green", lty = 4)
axis(side=1, at=seq(-1, 1, by = 0.5), labels = T)
axis(side=1, at=seq(-1, 1, by = 0.25), labels = F) 





plot(x, sapply(x, qs), 
     ylab = "", 
     xlab = "x", 
     main = "Quadratic Spectral", type = "l",
     xaxt="n",
     ylim = range_y)

lines(x, sapply(x, lugsail, 
                lugsail_parameters = list(r = 2, c = 1/2), 
                the_kernel= qs), 
      col = "red", lty = 2)
lugsail_parameters <- get_lugsail_parameters(200, 2, method = "Adaptive")
lines(x, sapply(x, lugsail, 
                lugsail_parameters =lugsail_parameters, 
                the_kernel = qs), 
      col = "blue", lty = 3)
lugsail_parameters <- get_lugsail_parameters(big_T = 200, q=2, method = "Over")
lines(x, sapply(x, lugsail, 
                lugsail_parameters =lugsail_parameters, 
                the_kernel = parzen), 
      col = "green", lty = 4)
axis(side=1, at=seq(-1, 1, by = 1), labels = T)
axis(side=1, at=seq(-1, 1, by = 0.25), labels = F)
dev.off()

