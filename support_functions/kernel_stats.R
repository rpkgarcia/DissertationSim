
all_c2 <- list(bartlett = list(mother = 2/3, zero = 1.333333), 
               parzen = list(mother = 151/280, zero = 0.6024786), 
               qs = list(mother = 1, zero = 1.288037))


all_c1 <- list(bartlett = list(mother = 1, zero = 3/2), 
               parzen = list(mother = 3/4, zero = 0.8007812), 
               qs = list(mother = 1.25, zero = 1.519851))

all_q <- list(bartlett= 1, 
              parzen = 2, 
              qs =2)

all_g_q <- list(bartlett= 1, 
                parzen = 6, 
                qs =1.42)


all_w_q <- function(q, rho){
  if(q ==1){
    w_q <- (1-rho^2)
  } else if(q == 2){
    w_q <- (1-rho)^2
  }
  
  w_q <- 2*rho/w_q
  return(w_q)
}
