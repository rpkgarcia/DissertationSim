
# Support Functions  ------------------------------------------------------

url <- "https://raw.githubusercontent.com/rpkgarcia/DissertationSim/main/support_functions/"

source(paste(url, "est_autocov.R", sep = ""))
source(paste(url, "kernels.R", sep = ""))
source(paste(url, "fitted_cv.R", sep = ""))

url <- "https://raw.githubusercontent.com/rpkgarcia/DissertationSim/main/simulations"
source(paste(url, "multi_sim_f_stat.R", sep = ""))
source(paste(url, "optimal_b_zerolug.R"))

# Data Set ----------------------------------------------------------------


# https://fred.stlouisfed.org/series/GDP#0
# GDP Billions of Dollars, Seasonally Adjusted Annual Rate

setwd(paste(getwd(), "/Applications", sep = ""))
GDP <- read.csv("GDP.csv")
big_T <- nrow(GDP)
GDP$log <- c(NA, log(GDP$GDP[2:big_T])- log(GDP$GDP[1:(big_T -1)]))
GDP <- GDP[-1, ]

plot(GDP$log, type = "l")


# Model  ------------------------------------------------------------------

# Raw growth rate 
# y_t = log(x_t) - log(x_{t-1})

# Model 
# y_t ~ theta + e_t 
# e_t <- rho*y_t + eps_t 
# eps_t is iid random error. 


# Test if growth rate is around 2%, or 0.02
# That is, we are testing if the intercept is 0.02
# We are NOT testing rho 
# Last 20 years, upto 2020Q3
#GDP <- GDP[(292-80):292, ]
plot(GDP$log, type = "l")
abline(h = 0.02, col = "red")

null_values <- 0.02
fit <- arima(GDP$log , order = c(1, 0, 0))
d <- 1
big_T <- length(u)

# Get residuals, center them  (centered kernel)
u <- fit$residuals
u <- u - mean(u)

# What value b should we use?
b_opt <- optimal_b_zero(rho=fit$coef[1], alpha=0.05, 
                        big_T= big_T, kernel=1, mother=2, d=1)
M <- b_opt*big_T


# ------- AutoCovariance Matrices  -------
# [#] the lag (0, ..., big_T-1)
all_autocovariances <-sapply(0:(big_T-1), R, the_sim_data = matrix(u))
all_autocovariances <- matrix(all_autocovariances)

# ------- Zero-Lugsail Estimator  -------
lug_para <- get_lugsail_parameters(big_T, q= 1, method = "Zero")

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
B <- floor(1/b)
scale <- (B-1:d + 1)/B

# Calculate F*-statistics for d = 1, ...., d for each simulation 
omega_hat_inv = chol2inv(omega_hat)
num <- (fit$coef[2]- null_values)
F_stat = scale*(num%*% omega_hat_inv %*%num)/d

omega_hat
sqrt(F_stat)
