# -------------------------------------------------------------------------
# Libraries and data sets
library(ggplot2)

# -------------------------------------------------------------------------
# Parameters 
d <- 5
z <- qchisq(.95, d)
try_T <- seq(100, 10000, by = 500)
names(try_T) <- try_T
try_rho <- seq(0.05, 0.95, by= 0.05)
names(try_rho) <- try_rho
c2 <- 1.33 # Bartlett lugsail 
c1 <- 1.4  # Bartlett lugsail 
delta_sq <- 5.67

# -------------------------------------------------------------------------
# Function that generates the best b from my analytical derivation. 

b_opt_func <- function(the_rho, the_big_T){
  num <- log(-c2*(d + z)*(1+ the_rho)/(4*the_rho^2*the_big_T*log(the_rho)))
  den <- the_big_T*log(the_rho)
  b <- num/den
  return(b)
}

b_opt <- sapply(try_rho, function(the_rho){
  sapply(try_T, function(the_big_T){
    b <- b_opt_func(the_rho, the_big_T)
    return(b)
  })
})

# -------------------------------------------------------------------------
# Function for optim 

prob_null <- function(the_b, the_rho, the_big_T, neg = -1){
  pdf_term <- dchisq(z, d, ncp = delta_sq)*z
  cdf_term <- pchisq(z, d, ncp = delta_sq)
  p <- c2*.5*(d + z)*the_b + 2*the_rho^2*the_rho^(the_b*the_big_T)/(1+ the_rho)
  p <- cdf_term - pdf_term*p
  return(p*neg)
}

# optim(par=0.05, fn = prob_null,method = "Brent", the_rho = 0.7, 
#       the_big_T = 200, lower = 0, upper = 1 )



b_optim <- sapply(try_rho, function(the_rho){
  sapply(try_T, function(the_big_T){
    b <- optim(par=0.05, fn = prob_null,method = "Brent", the_rho = the_rho, 
          the_big_T = the_big_T, lower = 0, upper = 1 )
    return(b$par)
  })
})

optim(par=0.05, fn = prob_null,method = "Brent", the_rho = 0.25, 
      the_big_T = 200, lower = 0, upper = 1 )$par
b_opt_func(0.25, 200)

# -------------------------------------------------------------------------
# Making Plots: See where my rull is negative  

map_me <- b_opt
# Plot data 
my_rows <- nrow(map_me)
my_cols <- ncol(map_me)
gg_me <- data.frame(map= c(map_me), 
                    map_TF = as.numeric(c(map_me)>0),  # 1 is good 
                    the_T = rep(try_T, times = length(try_rho)), 
                    the_rho = rep(try_rho, each = length(try_T)))

# Plot of Scenarios. 
ggplot(gg_me, aes(the_T, the_rho, fill = map_TF))+
  geom_tile(color = "black")+
  scale_fill_gradient2(low="red", high="white", mid = "pink", midpoint = 0) +
  labs(x = "Sample Size", y = "Correlation Coef.", 
       title = "Red = Bad, White = Good")

# Red indicates where bandwidth rule results in a negative bandwidth. 


# I tried this with both optim and my derivatives and we do have a problem. 
# We get a bandwidth impractically small, or zero. 


# -------------------------------------------------------------------------
# Making Plots: See where b_opt and b_optim differ by alot 

map_me <- b_opt- b_optim 
# Plot data 
my_rows <- nrow(map_me)
my_cols <- ncol(map_me)
gg_me <- data.frame(map= c(map_me), 
                    map_TF = as.numeric(c(map_me)>0.02),  # 1 is bad, big difference
                    the_T = rep(try_T, times = length(try_rho)), 
                    the_rho = rep(try_rho, each = length(try_T)))

# Plot of Scenarios. 
ggplot(gg_me, aes(the_T, the_rho, fill = map_TF))+
  geom_tile(color = "black")+
  scale_fill_gradient2(low="red", high="white", mid = "pink", midpoint = 0) +
  labs(x = "Sample Size", y = "Correlation Coef.", 
       title = "Red = Bad, White = Good")

# Red indicates where bandwidth rule results in a negative bandwidth. 


# I tried this with both optim and my derivatives and we do have a problem. 
# We get a bandwidth impractically small, or zero. 


# -------------------------------------------------------------------------
# The next big thing to try is sun's representation for the null prob and
# compare it to my representation to the null probability. 


sun_prob_null <- function(the_b, the_rho, the_big_T, neg = -1){
  pdf_term <- the_b*c2*dchisq(z, (d+2), ncp = delta_sq)*z*delta_sq/2
  cdf_term <- pchisq(z, d, ncp = delta_sq)
  #p <- c2*.5*(d + z)*the_b + 2*the_rho^2*the_rho^(the_b*the_big_T)/(1+ the_rho)
  p <- cdf_term + pdf_term
  return(p*neg)
}

# Pick some arbitrary b value that is reasonable 
some_b <- 0.07 

# See the probability under null that I have 
null_probs <- sapply(try_rho, function(the_rho){
  sapply(try_T, function(the_big_T){
    prob <- prob_null(the_b  = some_b , the_rho = the_rho, the_big_T = the_big_T, 
                      neg = 1)
    return(prob)
  })
})


# See the probability under null that Sun has  
sun_null_probs <- sapply(try_rho, function(the_rho){
  sapply(try_T, function(the_big_T){
    prob <- sun_prob_null(the_b  = some_b , the_rho = the_rho,
                          the_big_T = the_big_T, neg = 1)
    return(prob)
  })
})
# -------------------------------------------------------------------------
# Plot the difference between sun's null prob and mine. 

diff_null_probs <- null_probs - sun_null_probs

map_me <- diff_null_probs
# Plot data 
my_rows <- nrow(map_me)
my_cols <- ncol(map_me)
gg_me <- data.frame(map= c(map_me), 
                    map_TF = as.numeric(c(map_me)>0),  # 1 is good 
                    the_T = rep(try_T, times = length(try_rho)), 
                    the_rho = rep(try_rho, each = length(try_T)))


# Plot of Scenarios. 
ggplot(gg_me, aes(the_T, the_rho, fill = map))+
  geom_tile(color = "black")+
  scale_fill_gradient2(low="blue", high="red", mid = "white", midpoint = 0) +
  labs(x = "Sample Size", y = "Correlation Coef.", 
       title = "Red = Bad, White = Good")
# -------------------------------------------------------------------------
# Observe how the null prob changes for each rho 
# 
# for(the_rho in 1:length(try_rho)){
#   plot( try_T, sun_null_probs[,the_rho ],type = "l", 
#        main = paste("rho =",try_rho[the_rho]), 
#        ylim = c(-0.35, .5), 
#        ylab = "Type 2 Error", 
#        xlab = "Sample Size (T)")
#   lines( try_T,null_probs[,the_rho], col = "red")
#   Sys.sleep(1)
# }

# -------------------------------------------------------------------------
# Observe how the null prob changes for b, for a given T and Rho  


# Pick some arbitrary b value that is reasonable 
try_b <- seq(0.001, 0.12, by = 0.001)

# See the probability under null that I have 
null_probs <- sapply(try_b, function(some_b){
    results <- sapply(try_rho, function(the_rho){
      sapply(try_T, function(the_big_T){
        prob <- prob_null(the_b  = some_b , the_rho = the_rho, the_big_T = the_big_T, 
                          neg = 1)
        return(prob)
      })
    })
    return(results)
})

null_probs <- data.frame(null_probs)
colnames(null_probs) <- try_b 
null_probs$rho <-rep(try_rho, each = length(try_T))
null_probs$big_T <- rep(try_T, times = length(try_rho))

# See the probability under null that Sun has  
sun_null_probs  <- sapply(try_b, function(some_b){
  results <- sapply(try_rho, function(the_rho){
    sapply(try_T, function(the_big_T){
      prob <- sun_prob_null(the_b  = some_b , the_rho = the_rho, the_big_T = the_big_T, 
                        neg = 1)
      return(prob)
    })
  })
  return(results)
})

sun_null_probs <- data.frame(sun_null_probs)
colnames(sun_null_probs) <- try_b 
sun_null_probs$rho <-rep(try_rho, each = length(try_T))
sun_null_probs$big_T <- rep(try_T, times = length(try_rho))

some_plots <- sample(1:nrow(null_probs), 20)

# for(row_index in sort(some_plots)){
#   plot(try_b, sun_null_probs[row_index,1:length(try_b)],
#        type = "l", 
#        main = paste("rho =",sun_null_probs$rho[row_index],
#                     "big_T =", sun_null_probs$big_T[row_index]), 
#        ylim = c(-0.5, .5), 
#        ylab = "Type 2 Error", 
#        xlab = "Bandwidth (b)")
#   lines(try_b, null_probs[row_index,1:length(try_b)], col = "red")
#   sun_min_index <- which.min(sun_null_probs[row_index,1:length(try_b)])
#   min_index <- which.min(null_probs[row_index,1:length(try_b)])
#   points(c(try_b[sun_min_index], 
#            try_b[min_index]), 
#          c(sun_null_probs[row_index, sun_min_index],
#            null_probs[row_index, min_index]), 
#          col = c("black", "red")
#   )
#   
#   the_b <- optim(par=0.05, fn = prob_null,method = "Brent", 
#                  the_rho = sun_null_probs$rho[row_index], 
#                  the_big_T = sun_null_probs$big_T[row_index], 
#                  lower = 0, upper = 1 )
#   the_b_opt <- b_opt_func(sun_null_probs$rho[row_index], 
#                           sun_null_probs$big_T[row_index])
#   the_prob_b_opt <- prob_null(the_b_opt, 
#                               sun_null_probs$rho[row_index], 
#                               sun_null_probs$big_T[row_index], 1)
#   points(c(the_b$par, the_b_opt), 
#          c(-the_b$value,the_prob_b_opt), pch = c(8, 5), col = c("purple", "green"))
#   legend("bottomleft", 
#          c("My rule", "Optim", "Search"), 
#          col = c("purple", "green", "red"), 
#          pch = c(8, 5, 1))
#   Sys.sleep(1)
# }

# 
# row_index <- 65
# 
# plot(try_b, sun_null_probs[row_index,1:length(try_b)],
#      type = "l", 
#      main = paste("rho =",sun_null_probs$rho[row_index],
#                   "big_T =", sun_null_probs$big_T[row_index]), 
#      ylim = c(-0.35, .5), 
#      ylab = "Type 2 Error", 
#      xlab = "Bandwidth (b)")
# lines(try_b, null_probs[row_index,1:length(try_b)], col = "red")
# sun_min_index <- which.min(sun_null_probs[row_index,1:length(try_b)])
# min_index <- which.min(null_probs[row_index,1:length(try_b)])
# points(c(try_b[sun_min_index], 
#          try_b[min_index]), 
#        c(sun_null_probs[row_index, sun_min_index],
#          null_probs[row_index, min_index]), 
#        col = c("black", "red")
# )



# alternative rule  -------------------------------------------------------
alt_prob_null <- function(the_b, the_rho, the_big_T, neg = -1){
  pdf_term1 <- the_b*c2*dchisq(z, (d+2), ncp = delta_sq)*z*delta_sq/2
  pdf_term2 <- the_b*c2*dchisq(z, d, ncp = delta_sq)*z*(d-1)
  cdf_term <- pchisq(z, d, ncp = delta_sq)
  pdf_term3 <- dchisq(z, d, ncp = delta_sq)*z*2*the_rho^2*the_rho^(the_b*the_big_T)/(1+ the_rho)
  p <- cdf_term + pdf_term1 - pdf_term2 - pdf_term3
  return(p*neg)
}

# See the probability under alternative expression  
alt_null_probs  <- sapply(try_b, function(some_b){
  results <- sapply(try_rho, function(the_rho){
    sapply(try_T, function(the_big_T){
      prob <- alt_prob_null(the_b  = some_b , the_rho = the_rho, the_big_T = the_big_T, 
                            neg = 1)
      return(prob)
    })
  })
  return(results)
})

alt_null_probs <- data.frame(alt_null_probs)
colnames(alt_null_probs) <- try_b 
alt_null_probs$rho <-rep(try_rho, each = length(try_T))
alt_null_probs$big_T <- rep(try_T, times = length(try_rho))


# for(row_index in sort(some_plots)){
#   plot(try_b, sun_null_probs[row_index,1:length(try_b)],
#        type = "l", 
#        main = paste("rho =",sun_null_probs$rho[row_index],
#                     "big_T =", sun_null_probs$big_T[row_index]), 
#        ylim = c(-0.5, .5), 
#        ylab = "Type 2 Error", 
#        xlab = "Bandwidth (b)")
#   lines(try_b, null_probs[row_index,1:length(try_b)], col = "red")
#   lines(try_b, alt_null_probs[row_index,1:length(try_b)], col = "orange", lty = 2)
#   sun_min_index <- which.min(sun_null_probs[row_index,1:length(try_b)])
#   min_index <- which.min(null_probs[row_index,1:length(try_b)])
#   points(c(try_b[sun_min_index], 
#            try_b[min_index]), 
#          c(sun_null_probs[row_index, sun_min_index],
#            null_probs[row_index, min_index]), 
#          col = c("black", "red")
#   )
#   
#   the_b <- optim(par=0.05, fn = prob_null,method = "Brent", 
#                  the_rho = sun_null_probs$rho[row_index], 
#                  the_big_T = sun_null_probs$big_T[row_index], 
#                  lower = 0, upper = 1 )
#   the_b_opt <- b_opt_func(sun_null_probs$rho[row_index], 
#                           sun_null_probs$big_T[row_index])
#   the_prob_b_opt <- prob_null(the_b_opt, 
#                               sun_null_probs$rho[row_index], 
#                               sun_null_probs$big_T[row_index], 1)
#   points(c(the_b$par, the_b_opt), 
#          c(-the_b$value,the_prob_b_opt), pch = c(8, 5), col = c("purple", "green"))
#   legend("bottomleft", 
#          c("My rule", "Optim", "Search", "alternate expression"), 
#          col = c("purple", "green", "red", "orange"), 
#          pch = c(8, 5, 1, NA), 
#          lty = c(NA, NA, NA, 2))
#   Sys.sleep(1)
# }



# -------------------------------------------------------------------------
# Add expression for Type 1 error rate

error1 <- function(the_b, the_rho, the_big_T){
  alpha <- 0.05
  pdf_term2 <- dchisq(z, d)*z*c1*the_b
  pdf_term1 <- dchisq(z, d)*z*2*the_rho^2*the_rho^(the_b*the_big_T)/(1+ the_rho)
  p <- alpha + pdf_term1 #- pdf_term2
  return(p)
}


error1_probs  <- sapply(try_b, function(some_b){
  results <- sapply(try_rho, function(the_rho){
    sapply(try_T, function(the_big_T){
      prob <- error1(the_b  = some_b , the_rho = the_rho, the_big_T = the_big_T)
      return(prob)
    })
  })
  return(results)
})

error1_probs <- data.frame(error1_probs)
colnames(error1_probs) <- try_b 
error1_probs$rho <-rep(try_rho, each = length(try_T))
error1_probs$big_T <- rep(try_T, times = length(try_rho))


generate_info <- function(row_index, return_values =F){
  plot(try_b, sun_null_probs[row_index,1:length(try_b)],
       type = "l", 
       main = paste("rho =",sun_null_probs$rho[row_index],
                    "big_T =", sun_null_probs$big_T[row_index]), 
       ylim = c(-0.5, .5), 
       ylab = "Type 2 Error", 
       xlab = "Bandwidth (b)")
  lines(try_b, null_probs[row_index,1:length(try_b)], col = "red")
  lines(try_b, alt_null_probs[row_index,1:length(try_b)], col = "orange", lty = 2)
  sun_min_index <- which.min(sun_null_probs[row_index,1:length(try_b)])
  min_index <- which.min(null_probs[row_index,1:length(try_b)])
  points(c(try_b[sun_min_index], 
           try_b[min_index]), 
         c(sun_null_probs[row_index, sun_min_index],
           null_probs[row_index, min_index]), 
         col = c("black", "red")
  )
  
  the_b <- optim(par=0.05, fn = prob_null,method = "Brent", 
                 the_rho = sun_null_probs$rho[row_index], 
                 the_big_T = sun_null_probs$big_T[row_index], 
                 lower = 0, upper = 1 )
  the_b_opt <- b_opt_func(sun_null_probs$rho[row_index], 
                          sun_null_probs$big_T[row_index])
  the_prob_b_opt <- prob_null(the_b_opt, 
                              sun_null_probs$rho[row_index], 
                              sun_null_probs$big_T[row_index], 1)
  points(c(the_b$par, the_b_opt), 
         c(-the_b$value,the_prob_b_opt), pch = c(5, 8), 
         col = c("green", "purple"))
  legend("bottomleft", 
         c("My rule", "Optim", "Search", "alternate expression"), 
         col = c("purple", "green", "red", "orange"), 
         pch = c(8, 5, 1, NA), 
         lty = c(NA, NA, NA, 2))
  
  # Type 1 Error plot 
  plot(try_b, error1_probs[row_index,1:length(try_b)], 
       main = "Error 1", 
       ylab = "Type 1 Error", 
       xlab = "Bandwidth (b)", 
       type = "l")
  abline(h = 0.05, lty = 2)
  points(c(try_b[sun_min_index], 
           try_b[min_index]), 
         c(error1_probs[row_index, sun_min_index],
           error1_probs[row_index, min_index]), 
         col = c("black", "red")
  )
  
  e1_b_opt <- error1(the_b_opt,
                     sun_null_probs$rho[row_index], 
                     sun_null_probs$big_T[row_index])
  e1_the_b <- error1(the_b$par,
                     sun_null_probs$rho[row_index],
                     sun_null_probs$big_T[row_index])
  
  points(c(the_b$par, the_b_opt), 
         c(e1_the_b, e1_b_opt), 
         pch = c(8, 5), col = c("purple", "green"))
  
  legend("topright", 
         c("My rule", 
           "Optim", 
           "Search", 
           "alternate expression"), 
         col = c("purple", "green", "red", "orange"), 
         pch = c(8, 5, 1, NA), 
         lty = c(NA, NA, NA, 2))
  
  type1_error <- c("Suns_power_min" = error1_probs[row_index, sun_min_index],
                   "Search_power_min" = error1_probs[row_index, min_index], 
                   "Rule_power_min" = e1_the_b, 
                   "Optim_power_min" = e1_b_opt)
  type2_error <- c(sun_null_probs[row_index, sun_min_index],
                   null_probs[row_index, min_index], 
                   the_prob_b_opt, -the_b$value)
  b_value <- c(try_b[sun_min_index], 
               try_b[min_index],
               the_b_opt, 
               the_b$par)
  type <- c("Suns_power_min",
            "Search_power_min",
            "My_rule", 
            "Optim_rule")
  results <- data.frame(type1_error,type2_error, b_value,type)
  if(return_values == T){
    print(results) 
  }
}

par(mfrow = c(1, 2))
some_plots <- sample(1:nrow(null_probs), 5)
for(row_index in sort(some_plots)){
  generate_info(row_index, T)
  Sys.sleep(1)
}

p<- 0.9
generate_info(341, T)
-c2*(d + z)*z/(2*log(p)*0.5)

par(mfrow = c(1, 1))


# New b rule based on Type 1 error  ---------------------------------------
the_tau <- 1.0000
the_tau <- 1+dchisq(z, d)/2
b_opt_error1 <- function(the_rho, the_big_T, tau = the_tau){
  num <- log(0.05*(tau -1)*(1+ the_rho)/( dchisq(z, d)*2*the_rho^2*z))
  den <- the_big_T*log(the_rho)
  b <- num/den
  return(b)
}



  b_opt_e1  <- sapply(try_rho, function(the_rho){
    sapply(try_T, function(the_big_T){
      prob <- b_opt_error1(the_rho = the_rho, 
                           the_big_T = the_big_T)
      return(prob)
    })
  })

# b_opt_e1 <- data.frame(b_opt_e1)
# colnames(b_opt_e1) <- try_b 
# b_opt_e1$rho <-rep(try_rho, each = length(try_T))
# b_opt_e1$big_T <- rep(try_T, times = length(try_rho))


generate_info2 <- function(row_index, return_values =F){
  plot(try_b, sun_null_probs[row_index,1:length(try_b)],
       type = "l", 
       main = paste("rho =",sun_null_probs$rho[row_index],
                    "big_T =", sun_null_probs$big_T[row_index]), 
       ylim = c(-1, 1), 
       ylab = "Type 2 Error", 
       xlab = "Bandwidth (b)")
  lines(try_b, null_probs[row_index,1:length(try_b)], col = "red")
  lines(try_b, alt_null_probs[row_index,1:length(try_b)], col = "orange", lty = 2)
  sun_min_index <- which.min(sun_null_probs[row_index,1:length(try_b)])
  min_index <- which.min(null_probs[row_index,1:length(try_b)])
  points(c(try_b[sun_min_index], 
           try_b[min_index]), 
         c(sun_null_probs[row_index, sun_min_index],
           null_probs[row_index, min_index]), 
         col = c("black", "red")
  )
  
  the_b <- optim(par=0.05, fn = prob_null,method = "Brent", 
                 the_rho = sun_null_probs$rho[row_index], 
                 the_big_T = sun_null_probs$big_T[row_index], 
                 lower = 0, upper = 1 )
  the_b_opt <- b_opt_func(sun_null_probs$rho[row_index], 
                          sun_null_probs$big_T[row_index])
  the_prob_b_opt <- prob_null(the_b_opt, 
                              sun_null_probs$rho[row_index], 
                              sun_null_probs$big_T[row_index], 1)
  points(c(the_b$par, the_b_opt), 
         c(-the_b$value,the_prob_b_opt), pch = c(5, 8), 
         col = c("green", "purple"))
  
  the_b_opt_e1 <- b_opt_error1(sun_null_probs$rho[row_index], 
                               sun_null_probs$big_T[row_index])
  the_prob_b_opt_e1 <- prob_null(the_b_opt_e1, 
                              sun_null_probs$rho[row_index], 
                              sun_null_probs$big_T[row_index], 1)
  points(c(the_b_opt_e1), the_prob_b_opt_e1, 
         pch = 10, col = "blue")
  
  legend("bottomleft", 
         c("My rule", "Optim", "Search", "alternate expression", 
           "Error1-Optimal"), 
         col = c("purple", "green", "red", "orange", "blue"), 
         pch = c(8, 5, 1, NA, 10), 
         lty = c(NA, NA, NA, 2))
  
  # Type 1 Error plot 
  plot(try_b, error1_probs[row_index,1:length(try_b)], 
       main = "Error 1", 
       ylab = "Type 1 Error", 
       xlab = "Bandwidth (b)", 
       type = "l")
  abline(h = 0.05, lty = 2)
  abline(h = 0.0298+0.05, lty = 3, col = "grey")
  points(c(try_b[sun_min_index], 
           try_b[min_index]), 
         c(error1_probs[row_index, sun_min_index],
           error1_probs[row_index, min_index]), 
         col = c("black", "red")
  )
  
  e1_b_opt <- error1(the_b_opt,
                     sun_null_probs$rho[row_index], 
                     sun_null_probs$big_T[row_index])
  e1_the_b <- error1(the_b$par,
                     sun_null_probs$rho[row_index],
                     sun_null_probs$big_T[row_index])
  e1_b_opt_e1 <- error1(the_b_opt_e1,
                      sun_null_probs$rho[row_index],
                        sun_null_probs$big_T[row_index])
  
  points(c(the_b$par, the_b_opt, the_b_opt_e1), 
         c(e1_the_b, e1_b_opt, e1_b_opt_e1), 
         pch = c(8, 5, 10), col = c("purple", "green", "blue"))
  
  legend("topright", 
         legend = c("My rule", 
           "Optim", 
           "Search", 
            "error 1"), 
         col = c("purple", "green", "red",  "blue"), 
         pch = c(8, 5, 1,10))
  
  type1_error <- c("Suns_power_min" = error1_probs[row_index, sun_min_index],
                   "Search_power_min" = error1_probs[row_index, min_index], 
                   "Rule_power_min" = e1_the_b, 
                   "Optim_power_min" = e1_b_opt, 
                   "Error_1_rule"= e1_b_opt_e1)
  type2_error <- c(sun_null_probs[row_index, sun_min_index],
                   null_probs[row_index, min_index], 
                   the_prob_b_opt, -the_b$value, 
                   the_prob_b_opt_e1)
  b_value <- c(try_b[sun_min_index], 
               try_b[min_index],
               the_b_opt, 
               the_b$par, 
               the_b_opt_e1)
  type <- c("Suns_power_min",
            "Search_power_min",
            "My_rule", 
            "Optim_rule", 
            "Error_1_rule")
  results <- data.frame(type1_error,type2_error, b_value,type)
  if(return_values == T){
    print(results) 
    print(c("rho" = sun_null_probs$rho[row_index],
          "big_T" = sun_null_probs$big_T[row_index]))
  }
}

par(mfrow = c(1, 2))
some_plots <- sample(1:nrow(null_probs), 10)
for(row_index in sort(some_plots)){
  generate_info2(row_index, T)
  Sys.sleep(1)
}
par(mfrow = c(1, 1))



# -------------------------------------------------------------------------
c2_mother <- 2/3
g_q_mother <- 1

# ADD Sun expressions using standard Bartlett 
sun_prob_null2 <- function(the_b, the_rho, the_big_T){
  w_q <- 2*the_rho/(1- the_rho^2)
  pdf_term <- the_b*c2*dchisq(z, (d+2), ncp = delta_sq)*z*delta_sq/2
  cdf_term <- pchisq(z, d, ncp = delta_sq)
  pdf_term2 <- g_q_mother*w_q*dchisq(z, d, ncp = delta_sq)*z*(the_b*the_big_T)^(-1)
  p <- cdf_term + pdf_term - pdf_term2
  return(p)
}

# See the probability under alternative expression  
sun_prob_null2_values   <- sapply(try_b, function(some_b){
  results <- sapply(try_rho, function(the_rho){
    sapply(try_T, function(the_big_T){
      prob <- sun_prob_null2(the_b  = some_b , the_rho = the_rho, the_big_T = the_big_T)
      return(prob)
    })
  })
  return(results)
})

sun_prob_null2_values <- data.frame(sun_prob_null2_values)
colnames(sun_prob_null2_values) <- try_b 
sun_prob_null2_values$rho <-rep(try_rho, each = length(try_T))
sun_prob_null2_values$big_T <- rep(try_T, times = length(try_rho))


# Suns expresssion for type 1 
error1_sun <- function(the_b, the_rho, the_big_T){
  alpha <- 0.05
  w_q <- 2*the_rho/(1- the_rho^2)
  pdf_term3 <- g_q_mother*w_q*dchisq(z, d, ncp = delta_sq)*z*(the_b*the_big_T)^(-1)
  p <- alpha + pdf_term3
  return(p)
}


error1_sun_values  <- sapply(try_b, function(some_b){
  results <- sapply(try_rho, function(the_rho){
    sapply(try_T, function(the_big_T){
      prob <- error1_sun(the_b  = some_b , the_rho = the_rho, the_big_T = the_big_T)
      return(prob)
    })
  })
  return(results)
})


error1_sun_values <- data.frame(error1_sun_values)
colnames(error1_sun_values) <- try_b 
error1_sun_values$rho <-rep(try_rho, each = length(try_T))
error1_sun_values$big_T <- rep(try_T, times = length(try_rho))


# Get Sun's b optimal rule 
sun_b_rule <- function(the_rho, the_big_T, tau = the_tau){
  w_q <- 2*the_rho/(1- the_rho^2)
  num <- g_q_mother*w_q*dchisq(z, d)*z
  den <- (tau - 1)*0.05
  b <- (num/den)/the_big_T
  return(b)
}


sun_b_bartlett <- sapply(try_rho, function(the_rho){
  sapply(try_T, function(the_big_T){
    b <- sun_b_rule(the_rho, the_big_T)
    return(b)
  })
})



# Plot it all 



generate_info3 <- function(row_index, return_values =F){
  the_rho <- sun_null_probs$rho[row_index]
  the_big_T <-sun_null_probs$big_T[row_index]
  
  plot(try_b, sun_null_probs[row_index,1:length(try_b)],
       type = "l", 
       main = paste("rho =",sun_null_probs$rho[row_index],
                    "big_T =", sun_null_probs$big_T[row_index]), 
       ylim = c(-0.5, .5), 
       ylab = "Type 2 Error", 
       xlab = "Bandwidth (b)")
  lines(try_b, null_probs[row_index,1:length(try_b)], col = "red")
  lines(try_b, alt_null_probs[row_index,1:length(try_b)], col = "orange", lty = 2)
  sun_min_index <- which.min(sun_null_probs[row_index,1:length(try_b)])
  min_index <- which.min(null_probs[row_index,1:length(try_b)])
  points(c(try_b[sun_min_index], 
           try_b[min_index]), 
         c(sun_null_probs[row_index, sun_min_index],
           null_probs[row_index, min_index]), 
         col = c("black", "red")
  )
  
  the_b <- optim(par=0.05, fn = prob_null,method = "Brent", 
                 the_rho = sun_null_probs$rho[row_index], 
                 the_big_T = sun_null_probs$big_T[row_index], 
                 lower = 0, upper = 1 )
  the_b_opt <- b_opt_func(sun_null_probs$rho[row_index], 
                          sun_null_probs$big_T[row_index])
  the_prob_b_opt <- prob_null(the_b_opt, 
                              sun_null_probs$rho[row_index], 
                              sun_null_probs$big_T[row_index], 1)
  points(c(the_b$par, the_b_opt), 
         c(-the_b$value,the_prob_b_opt), pch = c(5, 8), 
         col = c("green", "purple"))
  
  the_b_opt_e1 <- b_opt_error1(sun_null_probs$rho[row_index], 
                               sun_null_probs$big_T[row_index])
  the_prob_b_opt_e1 <- alt_prob_null(the_b_opt_e1, 
                                 sun_null_probs$rho[row_index], 
                                 sun_null_probs$big_T[row_index], 1)
  points(c(the_b_opt_e1), the_prob_b_opt_e1, 
         pch = 10, col = "blue")
  
  
  # Add Sun - Mother Stuff 
  lines(try_b, sun_prob_null2_values[row_index,1:length(try_b)], col = "pink")
  sun_b_mother <- sun_b_rule(the_rho, the_big_T)
  sun_b_mother_prob <- sun_prob_null2(sun_b_mother, the_rho, the_big_T)
  points(sun_b_mother, sun_b_mother_prob, 
         pch = 12, col = "pink")
  
  legend("bottomright", 
         c("My rule", "Optim", "Search", "alternate expression", 
           "Error1-Optimal", "Sun-Mother"), 
         col = c("purple", "green", "red", "orange", "blue", "pink"), 
         pch = c(8, 5, 1, NA, 10, 12), 
         lty = c(NA, NA, NA, 2, NA, 1))
  
  # Type 1 Error plot 
  plot(try_b, error1_probs[row_index,1:length(try_b)], 
       main = "Error 1", 
       ylab = "Type 1 Error", 
       xlab = "Bandwidth (b)", 
       type = "l", 
       ylim = c(0.05, 0.062))
  abline(h = 0.05, lty = 2)
  abline(h = 0.0298+0.05, lty = 3, col = "grey")
  points(c(try_b[sun_min_index], 
           try_b[min_index]), 
         c(error1_probs[row_index, sun_min_index],
           error1_probs[row_index, min_index]), 
         col = c("black", "red")
  )
  
  e1_b_opt <- error1(the_b_opt,
                     sun_null_probs$rho[row_index], 
                     sun_null_probs$big_T[row_index])
  e1_the_b <- error1(the_b$par,
                     sun_null_probs$rho[row_index],
                     sun_null_probs$big_T[row_index])
  e1_b_opt_e1 <- error1(the_b_opt_e1,
                        sun_null_probs$rho[row_index],
                        sun_null_probs$big_T[row_index])
  
  points(c(the_b$par, the_b_opt, the_b_opt_e1), 
         c(e1_the_b, e1_b_opt, e1_b_opt_e1), 
         pch = c(8, 5, 10), col = c("purple", "green", "blue"))
  
  
  ### Add sun Stuff 
  lines(try_b, error1_sun_values[row_index,1:length(try_b)], 
        col = "pink", lty = 2)
  e1_sun <- error1_sun(sun_b_mother, 
                       sun_null_probs$rho[row_index], 
                       sun_null_probs$big_T[row_index])
  points(sun_b_mother, e1_sun, pch = 12, col = "pink", cex = 4)
  
  legend("topright", 
         legend = c("My rule", 
                    "Optim", 
                    "Search", 
                    "error 1"), 
         col = c("purple", "green", "red",  "blue", "pink"), 
         pch = c(8, 5, 1,10, 12))
  
  type1_error <- c("Suns_power_min" = error1_probs[row_index, sun_min_index],
                   "Search_power_min" = error1_probs[row_index, min_index], 
                   "Rule_power_min" = e1_the_b, 
                   "Optim_power_min" = e1_b_opt, 
                   "Error_1_rule"= e1_b_opt_e1, 
                   "Sun-Mother" = e1_sun)
  type2_error <- c(sun_null_probs[row_index, sun_min_index],
                   null_probs[row_index, min_index], 
                   the_prob_b_opt, -the_b$value, 
                   the_prob_b_opt_e1, 
                   sun_b_mother_prob)
  b_value <- c(try_b[sun_min_index], 
               try_b[min_index],
               the_b_opt, 
               the_b$par, 
               the_b_opt_e1, 
               sun_b_mother)
  type <- c("Suns_power_min",
            "Search_power_min",
            "My_rule", 
            "Optim_rule", 
            "Error_1_rule", 
            "Sun_Mother")
  results <- data.frame(type1_error,type2_error, b_value,type)
  if(return_values == T){
    print(results) 
    print(c("rho" = sun_null_probs$rho[row_index],
            "big_T" = sun_null_probs$big_T[row_index]))
  }
}


par(mfrow = c(1, 2))
some_plots <- sample(1:nrow(null_probs), 10)
for(row_index in sort(some_plots)){
  generate_info3(row_index, T)
  Sys.sleep(1)
}
par(mfrow = c(1, 1))


