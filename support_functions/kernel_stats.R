
all_c2 <- list(bartlett = list(mother = 2/3, zero = 1.333333), 
               parzen = list(mother = 151/280, zero = 0.6024786), 
               qs = list(mother = 1, zero = 1.288037))


all_c1 <- list(bartlett = list(mother = 1, zero = 3/2), 
               parzen = list(mother = 3/4, zero = 0.8007812), 
               qs = list(mother = 1.25, zero = 1.519851))

all_q <- list(bartlett= 1, 
              parzen = 2, 
              qs =2)

all_g_q <- list(bartlett= list(mother = 1, 
                               zero = 0, 
                               adaptive = NA, 
                               over = NA), 
                parzen = list(mother = 6, 
                              zero = 0, 
                              adaptive = NA, 
                              over = NA), 
                qs =list(mother = 1.42, 
                         zero = 0, 
                         adaptive = NA, 
                         over = NA))

para <- get_lugsail_parameters(200, 1, "Adaptive")
all_g_q$bartlett$adaptive <- (1-para$c*para$r)*all_g_q$bartlett$mother/(1-para$c)
para <- get_lugsail_parameters(200, 1, "Over")
all_g_q$bartlett$over <- (1-para$c*para$r)*all_g_q$bartlett$mother/(1-para$c)


para <- get_lugsail_parameters(200, 2, "Adaptive")
all_g_q$parzen$adaptive <- (1-para$c*para$r^2)*all_g_q$parzen$mother/(1-para$c)
para <- get_lugsail_parameters(200, 2, "Over")
all_g_q$parzen$over <- (1-para$c*para$r^2)*all_g_q$parzen$mother/(1-para$c)


para <- get_lugsail_parameters(200, 2, "Adaptive")
all_g_q$qs$adaptive <- (1-para$c*para$r^2)*all_g_q$qs$mother/(1-para$c)
para <- get_lugsail_parameters(200, 2, "Over")
all_g_q$qs$over <- (1-para$c*para$r^2)*all_g_q$qs$mother/(1-para$c)


all_w_q <- function(q, rho){
  if(q ==1){
    w_q <- (1-rho^2)
  } else if(q == 2){
    w_q <- (1-rho)^2
  }
  
  w_q <- 2*rho/w_q
  return(w_q)
}
