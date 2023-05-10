# Support Functions  ------------------------------------------------------

url <- "https://raw.githubusercontent.com/rpkgarcia/DissertationSim/main/support_functions/"

source(paste(url, "est_autocov.R", sep = ""))
source(paste(url, "kernels.R", sep = ""))
source(paste(url, "fitted_cv.R", sep = ""))

url <- "https://raw.githubusercontent.com/rpkgarcia/DissertationSim/main/simulations"
source(paste(url, "multi_sim_f_stat.R", sep = ""))
source(paste(url, "optimal_b_zerolug.R"))



# Data Set ----------------------------------------------------------------


#https://fred.stlouisfed.org/series/UNRATE#0
# Percent,Seasonally Adjusted (Monthly)
# Aprill 11, 2023

setwd(paste(getwd(), "/Applications/unemployment", sep = ""))
unrate <- read.csv("unrate.csv")
plot_me <- unrate$UNRATE[-(868:871)]
null_values <- 5.5


index <- seq(1,length(plot_me), by = 1)
#index <- 1:nrow(sp500)
plot(plot_me[index], type = "l")

fit_AR_1 <- arima(plot_me[index] , order = c(1, 0,0))
fit_AR_1
rho <- fit_AR_1$coef[1]

# ARMA MODEL FITTING 
AR_order <- rep(0:4, each = 4)
MA_order <- rep(0:4, times = 4)
model_info <- sapply(1:20, function(i) arima(plot_me[index],
                                             order= c(AR_order[i], 0, MA_order[i]))$aic)
model_info <- data.frame(AIC= model_info, AR = AR_order, MA = MA_order)
model_info <- model_info[order(model_info$AIC),]
head(model_info)

fit <- arima(plot_me[index] , order = c(model_info$AR[1], 0,model_info$MA[1]))
plot(fit$residuals)
fit
d <- 1

# Get residuals, center them  (centered kernel)
u <- fit$residuals
u <- u - mean(u)
big_T <- length(u)



# Bartlett Zero-Lugsail  --------------------------------------------------

# What value b should we use?
b_opt <- optimal_b_zero(rho=rho, alpha=0.05, 
                        big_T= big_T, kernel=1, mother=2, d=1)
M <- b_opt*big_T


# ------- AutoCovariance Matrices  -------
# [#] the lag (0, ..., big_T-1)
all_autocovariances <-sapply(0:(big_T-1), R, the_sim_data = matrix(u))
all_autocovariances <- matrix(all_autocovariances)
lug_para <- get_lugsail_parameters(big_T, q =1, method = "Zero")

# Weights matrix 
W <- rep(0, nrow(all_autocovariances))
nonzero_weights <- sapply(0:M/M, lugsail, 
                          the_kernel = bartlett, 
                          lugsail_parameters = lug_para)

# Everything beyond M has a zero weight
W[1:length(nonzero_weights)] <- nonzero_weights

# Calculate omega_hat 
omega_hat <- matrix(colSums(W * all_autocovariances)/big_T, d, d)

# Scaled factor for F-statistic,
B <- floor(1/b_opt)
scale <- (B-1:d + 1)/B

if(d == 1){
  scale <- 1
}

# Calculate F*-statistics for d = 1, ...., d for each simulation 
omega_hat_inv = chol2inv(omega_hat)
num <- (fit$coef[2]- null_values)
F_stat = scale*(num%*% omega_hat_inv %*%num)/d

omega_hat
if(d == 1){
  t_stat <- sqrt(F_stat)
}
sqrt(F_stat)

# ------- Test  -------
CV <- get_cv(d = d, alpha=0.05,
             the_kernel= "Bartlett", 
             is_lugsail = T, new_b = b_opt)
F_stat > CV # If true we reject the NULL 

zero_bartlett_lugsail <- c(CV = CV, 
                           F_stat = F_stat, 
                           Omega_hat = omega_hat)

# Bartlett Mother   --------------------------------------------------

# What value b should we use?
b_opt <- optimal_b_zero(rho=rho, alpha=0.05, 
                        big_T= big_T, kernel=1, mother=1, d=1)
M <- b_opt*big_T


# ------- AutoCovariance Matrices  -------
# [#] the lag (0, ..., big_T-1)
all_autocovariances <-sapply(0:(big_T-1), R, the_sim_data = matrix(u))
all_autocovariances <- matrix(all_autocovariances)
lug_para <- get_lugsail_parameters(big_T, q =1, method = "Mother")

# Weights matrix 
W <- rep(0, nrow(all_autocovariances))
nonzero_weights <- sapply(0:M/M, lugsail, 
                          the_kernel = bartlett, 
                          lugsail_parameters = lug_para)

# Everything beyond M has a zero weight
W[1:length(nonzero_weights)] <- nonzero_weights

# Calculate omega_hat 
omega_hat <- matrix(colSums(W * all_autocovariances)/big_T, d, d)

# Scaled factor for F-statistic,
B <- floor(1/b_opt)
scale <- (B-1:d + 1)/B

if(d == 1){
  scale <- 1
}

# Calculate F*-statistics for d = 1, ...., d for each simulation 
omega_hat_inv = chol2inv(omega_hat)
num <- (fit$coef[2]- null_values)
F_stat = scale*(num%*% omega_hat_inv %*%num)/d

omega_hat
if(d == 1){
  t_stat <- sqrt(F_stat)
}
sqrt(F_stat)

# ------- Test  -------
CV <- get_cv(d = d, alpha=0.05,
             the_kernel= "Bartlett", 
             is_lugsail = F, new_b = b_opt)
F_stat > CV # If true we reject the NUL


mother_bartlett_lugsail <- c(CV = CV, 
                           F_stat = F_stat, 
                           Omega_hat = omega_hat)



# ALL ---------------------------------------------------------------------
round(data.frame(mother_bartlett_lugsail, zero_bartlett_lugsail),4)

