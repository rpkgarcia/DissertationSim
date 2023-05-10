

# Page 461 of Priestly 
# Page 665 of Sun2014Lets 
all_c2 <- list(bartlett = list(mother = 2/3, 
                               zero = NA, 
                               adapt = NA, 
                               over = NA), 
               parzen = list(mother = 151/280, 
                             zero = NA, 
                             adapt = NA, 
                             over = NA), 
               qs = list(mother = 1, 
                         zero = NA, 
                         adapt = NA, 
                         over= NA))


# -------------------------------------------------------------------------
# Calculate c2 ------------------------------------------------------------
# -------------------------------------------------------------------------


# Bartlett ZeroLugsail  ~~~~~~~~~~~~~~~~~~~~~~~~~
support = function(x, M = 1, r= 2, c = 1/r){
  y = (1- abs(x/M))/(1-c) - c*(1- abs(x*r/M))/(1-c)
  y
}

bartlett_lug_sq = function(x, M=1, r=2 , c = 1/r, power = 2){
  w = ifelse(abs(x)<M, 
             ifelse(abs(x)<M/r,  support(x, r= r, c=c), 
                    (1-abs(x/M))/(1-c)), 0)
  return(w^power)
}
all_c2$bartlett$zero = integrate(bartlett_lug_sq, -1, 1)$value  # 1.333
para <- get_lugsail_parameters(200, 1, "Adaptive")
all_c2$bartlett$adapt = integrate(bartlett_lug_sq, -1, 1, r = para$r, c = para$c)$value  #1.57
all_c2$bartlett$over = integrate(bartlett_lug_sq, -1, 1, r = 2, c =(2/3))$value  # 2.333333

# Parzen ZeroLugsail  ~~~~~~~~~~~~~~~~~~~~~~~~~
# For values between 0 and 0.25
parzen1 <- function(x, M = 1, r = 2, c = .25){
  p1 <- 1 - 6*(x/M)^2 + 6*abs(x/M)^3
  p2 <- 1 - 6*(x*r/M)^2 + 6*abs(x*r/M)^3
  k_x <- (p1 - c*p2)/(1-c)
  return(k_x)
}

# For values between 0.25 and 0.5
parzen2 <- function(x, M = 1, r = 2, c = .25){
  p1 <- 1 - 6*(x/M)^2 + 6*abs(x/M)^3
  p2 <- 2*(1-abs(x/M))^3
  k_x <- (p1 - c*p2)/(1-c)
  return(k_x)
}

# For values between 0.5 and 1 
parzen3 <- function(x, M = 1, r = 2, c = .25){
  k_x <- 2*(1-abs(x/M))^3/(1-c)
  return(k_x)
}

parzen_lug_sq <- function(x, M=1, r=2 , c = .25){
  w <- ifelse(abs(x)<M, 
             ifelse(abs(x)<M/2,  
                    ifelse(abs(x)<M/(r*2), 
                           parzen1(x, M, r, c), parzen2(x, M, r, c)), 
                    parzen3(x, M, r, c)), 
             0)
  return(w^2)
}

all_c2$parzen$zero = integrate(parzen_lug_sq, -1, 1)$value  # 0.6024786
para <- get_lugsail_parameters(200, 2, "Adaptive")
all_c2$parzen$adapt = integrate(parzen_lug_sq, -1, 1, r = para$r, c = para$c)$value  # 0.618
para <- get_lugsail_parameters(200, 2, "Over")
all_c2$parzen$over = integrate(parzen_lug_sq, -1, 1, r = para$r, c =para$c)$value  #0.569


p1 <-integrate(parzen1, 0, 0.25, r = 1, c = 0)$value
p2<-integrate(parzen2, 0.25, 0.5, r = 1, c = 0)$value
p3<-integrate(parzen3, 0.5, 1, r = 1, c = 0)$value
c1_parzen <- 2*(p1+p2+p3)

p1 <-integrate(parzen1, 0, 0.25)$value
p2<-integrate(parzen2, 0.25, 0.5)$value
p3<-integrate(parzen3, 0.5, 1)$value
c1_zero_parzen <- 2*(p1+p2+p3)


para <- get_lugsail_parameters(200, 2, "Adaptive")
p1 <-integrate(parzen1, 0, 0.25, r = para$r, c = para$c)$value
p2<-integrate(parzen2, 0.25, 0.5, r = para$r, c = para$c)$value
p3<-integrate(parzen3, 0.5, 1, r = para$r, c = para$c)$value
c1_adapt_parzen <- 2*(p1+p2+p3)

para <- get_lugsail_parameters(200, 2, "Over")
p1 <-integrate(parzen1, 0, 0.25, r = para$r, c = para$c)$value
p2<-integrate(parzen2, 0.25, 0.5, r = para$r, c = para$c)$value
p3<-integrate(parzen3, 0.5, 1, r = para$r, c = para$c)$value
c1_over_parzen <- 2*(p1+p2+p3)

# QS ZeroLugsail  ~~~~~~~~~~~~~~~~~~~~~~~~~
# Same as above, but doesn't have 0 
qs_nonzero <- function(x){
  p1 <- sin(6*pi*x/5)/(6*pi*x/5)
  p2 <- cos(6*pi*x/5)
  p3 <- 25/(12*pi^2*x^2)
  k_x <- p3*(p1-p2) 
  
  k_x <- ifelse(x == 0 , 1, k_x)
  return(k_x)
}

qs_lug_sq <- function(x, M=1, r=2 , c = 1/r^(2), power = 2){
  w <- (qs_nonzero(x) - c*qs_nonzero(x*r))/(1-c)
  return(w^power)
}
all_c2$qs$zero = 2*integrate(qs_lug_sq, 0, 1)$value  # 1.288037
para <- get_lugsail_parameters(200, 2, "Adaptive")
all_c2$qs$adapt = integrate(qs_lug_sq, -1, 1, r = para$r, c = para$c)$value  # 0.618
para <- get_lugsail_parameters(200, 2, "Over")
all_c2$qs$over = integrate(qs_lug_sq, -1, 1, r = para$r, c =para$c)$value  #0.569






# -------------------------------------------------------------------------
# Calculate c1 ------------------------------------------------------------
# -------------------------------------------------------------------------

all_c1 <- list(bartlett = list(mother = 1, 
                               zero = 3/2, 
                               adapt = NA), 
               parzen = list(mother = 3/4, 
                             zero = c1_zero_parzen, 
                             adapt = c1_adapt_parzen, 
                             over = c1_over_parzen), 
               qs = list(mother = 1.25, 
                         zero = 1.519851, 
                         adapt = NA, 
                         over = NA ))


all_c1$qs$zero = 2*integrate(qs_lug_sq, 0, 1, power = 1)$value  # 1.288037
para <- get_lugsail_parameters(200, 2, "Adaptive")
all_c1$qs$adapt = integrate(qs_lug_sq, -1, 1, 
                            r = para$r, c = para$c, power =1 )$value  # 0.618
para <- get_lugsail_parameters(200, 2, "Over")
all_c1$qs$over = integrate(qs_lug_sq, -1, 1, 
                           r = para$r, c =para$c, power =1 )$value  #0.569


all_c1$bartlett$zero = integrate(bartlett_lug_sq, -1, 1, power = 1)$value  # 1.333
para <- get_lugsail_parameters(200, 1, "Adaptive")
all_c1$bartlett$adapt = integrate(bartlett_lug_sq, -1, 1,
                                  r = para$r, c = para$c, power = 1)$value  #1.57
all_c1$bartlett$over = integrate(bartlett_lug_sq, -1, 1, 
                                 r = 2, c =(2/3), power =1 )$value  # 2.333333


