# -------------------------------------------------------------------------
# Settings 
# -------------------------------------------------------------------------

input <- list(b = 0.07, the_kernel = "Bartlett")


try_b <-  seq(0.005, .99, by = 0.005)
index_b <- which(try_b == input$b)

# -------------------------------------------------------------------------
# Support Functions 
# -------------------------------------------------------------------------


# Get Critical Values 
fitted_model <- function(cv_matrix, alpha_level=0.05, m=1){
  chisq_cv <- qchisq(1-alpha_level, df = m)/m 
  try_b <- cv_matrix[,1]
  try_b <- as.numeric(gsub("b=", "", try_b))
  specific_cvs = cv_matrix[,m+1] - chisq_cv
  fit = lm(specific_cvs ~ 0+  poly(try_b, 3, raw = T))     
  return(fit)
}

# Get specific fitted/predicted value 
fitted_value <- function(fit, b, alpha = 0.05, m = 1){
  chisq_cv <- qchisq(1-alpha, df = m)/m
  fit_b <- sum(fit$coefficients*poly(b, 3, raw = T))
  fit_b <- chisq_cv + fit_b
  return(fit_b)
}


# -------------------------------------------------------------------------
# Distribution vs b Plot  
# -------------------------------------------------------------------------

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


# Bartlett Fits (need this for plot bounds)
file_name <- paste(CV_link, 
                   "Bartlett_", 
                   input$alpha, ".csv", sep = "")
bart_table <- read.csv(file_name)
bart_fit <- fitted_model(bart_table, alpha_level = alpha)
bart_cv <- fitted_value(bart_fit, as.numeric(input$b), alpha)

# Bartlett Lugsail Fits (need this for plot bounds)
file_name <- paste(CV_link, 
                   "Bartlett_Lugsail_", 
                   input$alpha, ".csv", sep = "")
bart_lug_table <- read.csv(file_name)
bart_lug_fit <- fitted_model(bart_lug_table, alpha_level = (alpha))
bart_lug_cv <- fitted_value(bart_lug_fit, as.numeric(input$b), alpha)


# Plot density curves
x_max <-max(10, bart_cv, bart_lug_cv)
plot(density(Mother[index,index_b], from = 0, to = 10), 
     xlim = xlim = c(0.01, x_max*1.03), ylim = c(0,1),
     main = "Distribution of Test Statistics", 
     ylab = the_ylab, xlab  = "F-Statistic")
mtext(cex=1, side = 3, "Given: kernel, d, bandwidth (b)")
curve(dchisq(x, 1), add = T, lty = 2, col = "grey", lwd = 2)
lines(density(Zero[index,index_b], from = 0, to = 10), 
      xlim = the_x_lim, col = "red", lty = 2)
lines(density(Adapt[index,index_b], from = 0, to = 10), 
      xlim = the_x_lim, col =  "blue", lty = 3)
lines(density(Over[index,index_b], from = 0, to = 10), 
      xlim = the_x_lim, col = "green", lty =4)




# Standard criticval value 
cv <- qchisq(1-alpha, 1)
lines(c(cv, cv), c(0, .20), col = "black", lty = 3, lwd = 4)
text(cv , .22, round(cv, 2))

# Bartlett CV
lines(c(bart_cv, bart_cv), c(0, .24), col = "red", lty = 3, lwd = 4)
text(bart_cv , .26, round(bart_cv, 2), col = "red")

# Bartlett_lugsail CV
lines(c(bart_lug_cv, bart_lug_cv), 
      c(0, .20), col = "blue", lty = 3, lwd = 4)
text(bart_lug_cv , .22, round(bart_lug_cv, 2), col = "blue")

# Legend 
legend("topright", 
       legend = c("Small-b", 
                  "Fixed-b: Bartlett", 
                  "Fixed-b: Lugsail"), 
       lty = c(2, 1, 1), 
       col = c("black", "red", "blue"), 
       lwd = 4)

# -------------------------------------------------------------------------
# CV for b Plot 
# -------------------------------------------------------------------------

# Alot of this came from "fitted_fixed_b.R"
alpha <- as.numeric(paste(".", input$alpha, sep = ""))
core_link <- "https://raw.githubusercontent.com/rpkgarcia/JobTalk/main/ShinyApps/Small_b_CV_distr/Bartlet_CV/"
bartlett_file <- paste(core_link, 
                       "Bartlett_", 
                       input$alpha, ".csv", sep = "")
bartlett_table <- read.csv(bartlett_file)
bartlett_lug_file <- paste(core_link, 
                           "Bartlett_Lugsail_", 
                           input$alpha, ".csv", sep = "")
bartlett_lug_table <- read.csv(bartlett_lug_file)


m <- 1
chisq_cv <- qchisq(1-alpha, df = m)/m 
try_b <- bartlett_table[,1]
try_b <- as.numeric(gsub("b=", "", try_b))
fit_bartlett <- fitted_model(bartlett_table, alpha_level = alpha)
fit_bartlett_lug <- fitted_model(bartlett_lug_table, alpha_level = alpha)

bartlett_fit <- fitted_value(fit_bartlett, 
                             as.numeric(input$b), alpha)
bartlett_lug_fit <- fitted_value(fit_bartlett_lug, 
                                 as.numeric(input$b), alpha) 
y_max <- max(10, bartlett_fit, bartlett_lug_fit)

plot(c(0, try_b), 
     c(chisq_cv, fit_bartlett$fitted.values + chisq_cv), 
     type = "l", 
     xlim = c(0, 0.5), ylim = c(0, y_max),
     col = "red", lwd = 4, xlab = "Bandwidth (b)", 
     ylab = "Critical Value", 
     main = c("Critical Values vs Bandwidth (b)"))
mtext( cex=1, side = 3,
       "Given: kernel, d, significance level (alpha)")
lines(c(0, try_b), 
      c(chisq_cv, fit_bartlett_lug$fitted.values + chisq_cv), 
      col = "blue", lwd = 4)

# My observed values 
abline(h = chisq_cv, lwd = 4, lty = 2)
abline(v = input$b, lwd = 4, lty = 3, col = "grey")

legend("bottomright", 
       legend = c("Small-b CV", 
                  "Fixed-b: Bartlett", 
                  "Fixed-b: Lugsail"), 
       lty = c(2, 1, 1, 3), 
       col = c("black", "red", "blue"), 
       lwd = 4)