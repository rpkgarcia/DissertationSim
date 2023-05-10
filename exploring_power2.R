# -------------------------------------------------------------------------
# Libraries and data sets
library(ggplot2)

# -------------------------------------------------------------------------
# Parameters 
d <-1
z <- qchisq(.95, d)
try_T <- seq(100, 1000, by = 100)
names(try_T) <- try_T
try_rho <- seq(0.05, 0.95, by= 0.05)
names(try_rho) <- try_rho
c2 <- 1.33 # Bartlett lugsail 
c1 <- 1.4  # Bartlett lugsail 
delta_sq <- 5.67
try_b <- seq(0.00001, .5, by = 0.01)


# -------------------------------------------------------------------------
# Sun's representation for Type 2 error 
# Zero Lugsail 

sun_error2_func <- function(the_b, the_rho, the_big_T){
  pdf_term <- the_b*c2*dchisq(z, (d+2), ncp = delta_sq)*z*delta_sq/2
  cdf_term <- pchisq(z, d, ncp = delta_sq)
  p <- cdf_term + pdf_term
  return(p)
}


# See the probability under alternative expression  
sun_error2<- sapply(try_b, function(some_b){
  results <- sapply(try_rho, function(the_rho){
    sapply(try_T, function(the_big_T){
      prob <- sun_error2_func(the_b  = some_b , the_rho = the_rho, the_big_T = the_big_T)
      return(prob)
    })
  })
  return(results)
})

sun_error2<- data.frame(sun_error2)
colnames(sun_error2) <- try_b 
sun_error2$rho <-rep(try_rho, each = length(try_T))
sun_error2$big_T <- rep(try_T, times = length(try_rho))




# -------------------------------------------------------------------------
# Sun's expression for Type 2 error 
# Mother kernel
c2_mother <- 2/3
g_q_mother <- 1

# ADD Sun expressions using standard Bartlett 
sun_error2_mother_func <- function(the_b, the_rho, the_big_T){
  w_q <- 2*the_rho/(1- the_rho^2)
  pdf_term <- the_b*c2_mother*dchisq(z, (d+2), ncp = delta_sq)*z*delta_sq/2
  cdf_term <- pchisq(z, d, ncp = delta_sq)
  pdf_term2 <- g_q_mother*w_q*dchisq(z, d, ncp = delta_sq)*z*(the_b*the_big_T)^(-1)
  p <- cdf_term + pdf_term - pdf_term2
  return(p)
}

# See the probability under alternative expression  
sun_error2_mother   <- sapply(try_b, function(some_b){
  results <- sapply(try_rho, function(the_rho){
    sapply(try_T, function(the_big_T){
      prob <- sun_error2_mother_func(the_b  = some_b , the_rho = the_rho, the_big_T = the_big_T)
      return(prob)
    })
  })
  return(results)
})

sun_error2_mother  <- data.frame(sun_error2_mother )
colnames(sun_error2_mother ) <- try_b 
sun_error2_mother$rho <-rep(try_rho, each = length(try_T))
sun_error2_mother$big_T <- rep(try_T, times = length(try_rho))

# -------------------------------------------------------------------------
# Sun's expression for Type 1 error 
# Mother kernel

# Suns expresssion for type 1 
sun_error1_mother_func <- function(the_b, the_rho, the_big_T){
  alpha <- 0.05
  w_q <- 2*the_rho/(1- the_rho^2)
  pdf_term3 <- g_q_mother*w_q*dchisq(z, d)*z*(the_b*the_big_T)^(-1)
  p <- alpha + pdf_term3
  return(p)
}


sun_error1_mother  <- sapply(try_b, function(some_b){
  results <- sapply(try_rho, function(the_rho){
    sapply(try_T, function(the_big_T){
      prob <- sun_error1_mother_func(the_b  = some_b , the_rho = the_rho, the_big_T = the_big_T)
      return(prob)
    })
  })
  return(results)
})


sun_error1_mother <- data.frame(sun_error1_mother)
colnames(sun_error1_mother) <- try_b 
sun_error1_mother$rho <-rep(try_rho, each = length(try_T))
sun_error1_mother$big_T <- rep(try_T, times = length(try_rho))

# -------------------------------------------------------------------------
# Sun's expression for Optimal b
# Mother kernel

the_tau <- 1 + .2

# Get Sun's b optimal rule 
sun_b_mother_func <- function(the_rho, the_big_T, tau = 1.2){
  w_q <- 2*the_rho/(1- the_rho^2)
  num <- g_q_mother*w_q*dchisq(z, d)*z
  den <- (tau - 1)*0.05
  b <- (num/den)/the_big_T
  return(b)
}


# sun_b_bartlett <- sapply(try_rho, function(the_rho){
#   sapply(try_T, function(the_big_T){
#     b <- sun_b_rule(the_rho, the_big_T)
#     return(b)
#   })
# })


# -------------------------------------------------------------------------
# My representation for Type 2 error 
# Zero Lugsail

error2_func <- function(the_b, the_rho, the_big_T){
  # pdf_term1 <- the_b*c2*dchisq(z, (d+2), ncp = delta_sq)*z*delta_sq/2
  # pdf_term2 <- the_b*c2*dchisq(z, d, ncp = delta_sq)*z*(d-1)
  # cdf_term <- pchisq(z, d, ncp = delta_sq)
  # pdf_term3 <- dchisq(z, d, ncp = delta_sq)*z*2*the_rho^2*the_rho^(the_b*the_big_T)/(1+ the_rho)
  # p <- cdf_term + pdf_term1 - pdf_term2 - pdf_term3
  
  
  cdf_term <- pchisq(z, d, delta_sq)
  pdf_term2 <- dchisq(z, d, delta_sq)*z*c1*the_b
  pdf_term1 <- dchisq(z, d, delta_sq)*z*2*the_rho^2*the_rho^(the_b*the_big_T)/(1+ the_rho)
  p <- cdf_term + pdf_term1 #- pdf_term2
  
  
  return(1-p)
}

# See the probability under alternative expression  
error2  <- sapply(try_b, function(some_b){
  results <- sapply(try_rho, function(the_rho){
    sapply(try_T, function(the_big_T){
      prob <- error2_func(the_b  = some_b , the_rho = the_rho, the_big_T = the_big_T)
      return(prob)
    })
  })
  return(results)
})

error2 <- data.frame(error2)
colnames(error2) <- try_b 
error2$rho <-rep(try_rho, each = length(try_T))
error2$big_T <- rep(try_T, times = length(try_rho))




# -------------------------------------------------------------------------
# My representation for Type 1 error 
# Zero Lugsail 

error1_func <- function(the_b, the_rho, the_big_T){
  alpha <- 0.05
  pdf_term2 <- dchisq(z, d)*z*c1*the_b
  pdf_term1 <- dchisq(z, d)*z*2*the_rho^2*the_rho^(the_b*the_big_T)/(1+ the_rho)
  p <- alpha + pdf_term1 #- pdf_term2
  return(p)
}


error1  <- sapply(try_b, function(some_b){
  results <- sapply(try_rho, function(the_rho){
    sapply(try_T, function(the_big_T){
      prob <- error1_func(the_b  = some_b , the_rho = the_rho, the_big_T = the_big_T)
      return(prob)
    })
  })
  return(results)
})

error1 <- data.frame(error1)
colnames(error1) <- try_b 
error1$rho <-rep(try_rho, each = length(try_T))
error1$big_T <- rep(try_T, times = length(try_rho))



# -------------------------------------------------------------------------
# The lower bound for b using my rule 
# Zero Lugsail 

the_tau <- 1.0000
the_tau <- 1+dchisq(z, d)/2

b_LB_func <- function(the_rho, the_big_T, tau = 1){
  tau <- 1 + 1/the_big_T
  num <- log(0.05*(1/the_big_T)*(1+ the_rho)/( dchisq(z, d)*2*the_rho^2*z))
  den <- the_big_T*log(the_rho)
  b <- num/den
  return(b)
}

# 
# 
# b_LB  <- sapply(try_rho, function(the_rho){
#   sapply(try_T, function(the_big_T){
#     prob <- b_LB_func(the_rho = the_rho, 
#                          the_big_T = the_big_T)
#     return(prob)
#   })
# })




# -------------------------------------------------------------------------
# Plot it all 



generate_info3 <- function(row_index, return_values =F){
  the_rho <- error1$rho[row_index]
  the_big_T <-error1$big_T[row_index]
  
  
  # ~~~~~  Type 1 Error ~~~~~  
  plot(try_b, sun_error1_mother[row_index,1:length(try_b)], 
       main = "Error 1", 
       ylab = "Type 1 Error", 
       xlab = "Bandwidth (b)", 
       type = "l", lty = 2, 
       col = "red")
  lines(try_b, error1[row_index,1:length(try_b)], 
        col = c("blue"), lty = 1)
  abline(h = 0.05, lty = 2)
  abline(h = 1/the_big_T +0.05, lty = 3, col = "grey")
  
  # Optimal b values 
  b_LB <- b_LB_func(the_rho, the_big_T)
  b_sun <- sun_b_mother_func(the_rho, the_big_T)
  
  # Error 1 Probs 
  e1_b_LB <- error1_func(b_LB, the_rho, the_big_T)
  e1_b_sun <- sun_error1_mother_func(b_sun, the_rho, the_big_T)
  
  # Plot points 
  points(c(b_LB, b_sun), 
         c(e1_b_LB, e1_b_sun), 
         col = c("blue", "red"), 
         pch = c(8, 10))
  
  legend("topright", 
         legend = c("Zero", 
                    "Mother (Sun)"), 
         col = c("blue", "red"), 
         pch = c(8,10), 
         lty = c(1, 2))
  
  # ~~~~~  Type 2 Error ~~~~~  
  plot(try_b, sun_error2[row_index,1:length(try_b)],
       type = "l", 
       main = paste("rho =",the_rho,
                    "big_T =", the_big_T), 
       ylim = c(0, 1),
       ylab = "Type 2 Error", 
       xlab = "Bandwidth (b)", lwd = 1, col = "red")
  lines(try_b, error2[row_index,1:length(try_b)], col = "blue", lwd =2 )
  lines(try_b, sun_error2_mother[row_index,1:length(try_b)], col = "red", lty = 2)
  abline(h = pchisq(z, 1, ncp = delta_sq))

  # Error 1 Probs 
  e2_b_LB <- error2_func(b_LB, the_rho, the_big_T)
  e2_b_sun <- sun_error2_mother_func(b_sun, the_rho, the_big_T)
  
  
  # Plot points 
  points(c(b_LB, b_sun), 
         c(e2_b_LB, e2_b_sun), 
         col = c("blue", "red"), 
         pch = c(8, 10))
  
  # legend("bottomright", 
  #        legend = c("Zero", 
  #                   "Mother (Sun)", 
  #                   "Zero (Sun)"), 
  #        col = c("blue", "red", "red"), 
  #        pch = c(8,10, NA), 
  #        lty = c(1, 2, 1))
  
  
  # ~~~~~  Return Values  ~~~~~  
  type1_error <- c("Zero" = e1_b_LB,
                   "Mother(Sun)" = e1_b_sun)
  type2_error <- c(e2_b_LB, 
                   e2_b_sun)
  b_value <- c(b_LB, 
               b_sun)
  results <- data.frame(type1_error,type2_error, b_value)
  if(return_values == T){
    print(results) 
    print(c("rho" = the_rho,
            "big_T" = the_big_T))
  }
}


par(mfrow = c(1, 2))
set.seed(62)
some_plots <- sample(1:nrow(error1), 10)
for(row_index in sort(some_plots)){
  generate_info3(row_index, T)
  Sys.sleep(1)
}
#par(mfrow = c(1, 1))



get_plot_specific<- function(the_rho, the_big_T){
  index <- which((round(error1$rho, 2)== the_rho) & (error1$big_T == the_big_T))
  generate_info3(index, T)
}

get_plot_specific(.7, 200)
get_plot_specific(.7, 100)
