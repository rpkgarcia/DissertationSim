QS_05 <- read.csv("~/Downloads/temp_CV/QS_05.csv")
QS_Lugsail_05 <- read.csv("~/Downloads/temp_CV/QS_Lugsail_05.csv")
Bartlett_05 <- read.csv("~/Downloads/temp_CV/Bartlett_05.csv")
Bartlett_Lugsail_05 <- read.csv("~/Downloads/temp_CV/Bartlett_Lugsail_05.csv")


adjust <- function(x){
  return((x))
}
index <- which(QS_05$X<0.4)

plot(QS_05$X[index], 
     adjust(QS_05$X1[index]), type = "l")
lines(QS_Lugsail_05$X[index], 
      adjust(QS_Lugsail_05$X1[index]), lty = 2)

lines(Bartlett_05$X[index], 
      adjust(Bartlett_05$X1[index]), lty = 1, col = "red")
lines(Bartlett_Lugsail_05$X[index], 
      adjust(Bartlett_Lugsail_05$X1[index]), lty = 2, col = "red")



# -------------------------------------------------------------------------

Bartlett_d1 <- read.csv("~/Downloads/rstudio-export (3)/Bartlett_d1.csv")
Bartlett_Lug_d1 <- read.csv("~/Downloads/rstudio-export (3)/Bartlett_Lugsail_d1.csv")
QS_d1 <- read.csv("~/Downloads/rstudio-export (3)/QS_d1.csv")
QS_Lug_d1 <- read.csv("~/Downloads/rstudio-export (3)/QS_Lugsail_d1.csv")

index <-which(Bartlett_d1[,15]<10)
index <- 1:nrow(Bartlett_d1)
#index <- sample(1:nrow(Bartlett_d1), 100000, replace = T)
plot(density(Bartlett_d1[index,15], from = 0, to = 10), col = "red")
#lines(density(Bartlett_d1[,15], from = 0, to = 10), lty =2, col = "red")
lines(density(Bartlett_Lug_d1[,15], from = 0, to = 10), lty =2, col = "red")
lines(density(QS_d1[,15], from = 0, to = 10), lty =1, col =1)
lines(density(QS_Lug_d1[,15], from = 0, to = 10), lty =2, col =1)
curve(dchisq(x, 1), add = T)

plot(density(distr_est_b0.07$Bartlett))

# -------------------------------------------------------------------------


Bartlett_05_big <- read.csv("~/Documents/GitHub/DissertationSim/fixed_b_asymptotics/Fixed_b_CV_Tables/Bartlett_05.csv")
Bartlett_Lug_05_big <- read.csv("~/Documents/GitHub/DissertationSim/fixed_b_asymptotics/Fixed_b_CV_Tables/Bartlett_Lugsail_05.csv")

plot(Bartlett_05_big$X[index], 
     adjust(Bartlett_05_big$X1[index]), type = "l")
lines(Bartlett_Lug_05_big$X[index], 
      adjust(Bartlett_Lug_05_big$X1[index]), lty = 2)

lines(Bartlett_05$X[index], 
      adjust(Bartlett_05$X1[index]), lty = 1, col = "red")
lines(Bartlett_Lugsail_05$X[index], 
      adjust(Bartlett_Lugsail_05$X1[index]), lty = 2, col = "red")



# -------------------------------------------------------------------------

m <- matrix(sample(1:9,9), nrow = 3)
m
solve(m)
chol2inv(m)
chol2inv(chol(m))
