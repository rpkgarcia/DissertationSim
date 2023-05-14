# -------------------------------------------------------------------------
# Settings 
# -------------------------------------------------------------------------

# Change this line 
input <- list(b = 0.10, the_kernel = "Bartlett", alpha = 0.10, d= 1)

# Do not change these lines 
try_b <-  seq(0.005, .99, by = 0.005)
index_b <- which(try_b == input$b)+ 1
#index_b <- 10  # for b = 0.05
setwd("~/Documents/GitHub/DissertationSim/fixed_b_asymptotics/1_dimensional")


# -------------------------------------------------------------------------
# Support Functions 
# -------------------------------------------------------------------------
# Read in all fitted values 
the_fits <- read.csv("Fixed_b_CV_Tables_New/fitted_CV.csv") 

# New_b is a vector of b values (could have length >1) of which 
# you want to find the CV for. 

get_cv <- function(d, alpha, the_kernel, lugsail, new_b){
  chisq_cv <-  qchisq(1-alpha, df = d)/d 
  
  # Pull out only the values you need 
  index <- the_fits$kernel == the_kernel & the_fits$lugsail == lugsail &
    the_fits$alpha == alpha & the_fits$d == d
  
  coefficients <- the_fits[index, c("beta1", "beta2","beta3")]
  intercept <-  the_fits[index, c("intercept")]
  
  # Fitted value of b 
  new_b<-data.frame(poly(new_b, 3, raw = T))
  cv_by_b <- apply(new_b, 1, function(x) sum(x*coefficients))
  cv_by_b <- cv_by_b + chisq_cv
  
  return(cv_by_b)
}


# -------------------------------------------------------------------------
# Save plot information ---------------------------------------------------
# -------------------------------------------------------------------------

start_plot <- function(plot_name){
  reso <- 600
  length <- 4
  png(paste(plot_name, ".png", sep = ""),
      units="in", 
      res = reso, 
      width = 8, 
      height = 4)
}


start_plot(paste(c("", "b", "alpha", "d"),
                c(input$the_kernel, input$b, input$alpha,input$d),
                sep = "", collapse = "_"))
par(mfrow = c(1, 2))
# -------------------------------------------------------------------------
# Distribution vs b Plot  
# -------------------------------------------------------------------------

# Obtaining Values for the distribution plots 
file_mother<- paste("Fixed_b_distribution_New/", input$the_kernel, "_", 
                   "Mother", ".csv", sep = "")
Mother <- read.csv(file_mother)

file_zero<- paste("Fixed_b_distribution_New/", input$the_kernel, "_", 
                    "Zero", ".csv", sep = "")
Zero <- read.csv(file_zero)

file_adapt<- paste("Fixed_b_distribution_New/", input$the_kernel, "_", 
                    "Adapt", ".csv", sep = "")
Adapt <- read.csv(file_adapt)

file_over<- paste("Fixed_b_distribution_New/", input$the_kernel, "_", 
                   "Over", ".csv", sep = "")
Over <- read.csv(file_over)


# Get CVs 
cv_mother <- get_cv(d = input$d, alpha = input$alpha,
                    the_kernel = input$the_kernel,
                    lugsail = "Mother", new_b = input$b)
cv_zero <- get_cv(d = input$d, alpha = input$alpha,
                    the_kernel = input$the_kernel,
                    lugsail = "Zero", new_b = input$b)
cv_adapt <- get_cv(d = input$d, alpha = input$alpha,
                    the_kernel = input$the_kernel,
                    lugsail = "Adapt", new_b = input$b)
cv_over <- get_cv(d = input$d, alpha = input$alpha,
                    the_kernel = input$the_kernel,
                    lugsail = "Over", new_b = input$b)

# -----   Plot Density Curves -----
# Plot density curves
x_max <-min(20, cv_over)
plot(density(Mother[index,index_b], from = 0, to = 20),
     xlim = c(0.01, x_max*1.03), ylim = c(0,1),
     main = "Density of Test Statistics",
     ylab="", xlab  = "F-Statistic")
mtext(cex=1, side = 3, "Given: kernel, d, bandwidth (b)")
curve(dchisq(x, 1), add = T, lty = 2, col = "darkgrey", lwd = 2)
lines(density(Zero[index,index_b], from = 0, to = 20),
      xlim = the_x_lim, col = "red", lty = 2)
lines(density(Adapt[index,index_b], from = 0, to = 20),
      xlim = the_x_lim, col =  "blue", lty = 3)
lines(density(Over[index,index_b], from = 0, to = 20),
      xlim = the_x_lim, col = "green", lty =4)




# -----   Add CV Lines -----
# Standard CV
cv_small_b <- qchisq(1-input$alpha, input$d)
lines(c(cv_small_b, cv_small_b), c(0, .40), col = "darkgrey", lty = 1, lwd = 1)
text(cv_small_b , .47, round(cv_small_b, 2), col = "black", cex = 0.75)

# Mother CV
lines(c(cv_mother, cv_mother), c(0, .20), col = "black", lty = 1, lwd = 1)
text(cv_mother , .27, round(cv_mother, 2), col = "black", cex =0.75)

# Zero CV
lines(c(cv_zero, cv_zero), c(0, .40), col = "red", lty = 1, lwd = 1)
text(cv_zero , .47, round(cv_zero, 2), col = "black", cex = 0.75)

# Adapt CV
lines(c(cv_adapt, cv_adapt),  c(0, .2), col = "blue", lty = 1, lwd = 1)
text(cv_adapt , .27, round(cv_adapt, 2), col = "black", cex = 0.75)

# Over CV
lines(c(cv_over, cv_over), c(0, .40), col = "green", lty = 1, lwd = 1)
text(cv_over, .47, round(cv_over, 2), col = "black", cex = 0.75)

# -----   Plot Legend -----
# legend("topright",
#        legend = c("Small-b",
#                   "Fixed-b: Bartlett",
#                   "Fixed-b: Lugsail"),
#        lty = c(2, 1, 1),
#        col = c("black", "red", "blue"),
#        lwd = 4)

# -------------------------------------------------------------------------
# CV for b Plot
# -------------------------------------------------------------------------
try_b_subset <- try_b[try_b <= 0.5]

fit_mother <- get_cv(d = input$d, alpha = input$alpha,
                     the_kernel = input$the_kernel,
                     lugsail = "Mother", new_b = try_b_subset)
fit_zero <- get_cv(d = input$d, alpha = input$alpha,
                     the_kernel = input$the_kernel,
                     lugsail = "Zero", new_b = try_b_subset)

fit_adapt <- get_cv(d = input$d, alpha = input$alpha,
                     the_kernel = input$the_kernel,
                     lugsail = "Adapt", new_b = try_b_subset)

fit_over <- get_cv(d = input$d, alpha = input$alpha,
                     the_kernel = input$the_kernel,
                     lugsail = "Over", new_b = try_b_subset)

y_max <- max(c(10, fit_over))

plot(c(0, 0.5),
     c(cv_small_b, cv_small_b),
     type = "l",
     xlim = c(0, 0.5), ylim = c(0, y_max*1.03), lty = 2,
     col = "darkgrey", lwd = 2, xlab = "Bandwidth (b)",
     ylab="",
     main = c("Critical Values"))
mtext(cex=1, side = 3, expression(paste("Given: kernel, d, ",alpha)))

lines(c(0, try_b_subset), c(cv_small_b, fit_mother), col = "black", lty = 1)
lines(c(0, try_b_subset), c(cv_small_b, fit_zero), col = "red", lty = 2)
lines(c(0, try_b_subset), c(cv_small_b, fit_adapt), col = "blue", lty = 3)
lines(c(0, try_b_subset), c(cv_small_b, fit_over), col = "green", lty = 4)

# My observed values
abline(v = input$b, lwd = 2, lty = 1, col = "black")

# legend("bottomright",
#        legend = c("Small-b CV",
#                   "Fixed-b: Bartlett",
#                   "Fixed-b: Lugsail"),
#        lty = c(2, 1, 1, 3),
#        col = c("black", "red", "blue"),
#        lwd = 4)

# Save plot with 1000 width, 400 hieght aspect ratio 
dev.off()

