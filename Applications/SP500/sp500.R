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

setwd(paste(getwd(), "/Applications/SP500", sep = ""))
sp500 <- read.csv("sp500.csv")
colnames(sp500)[ncol(sp500)] <- "Change"
sp500$Change <- as.numeric(gsub("%", "", sp500$Change))
sp500$Open <- as.numeric(gsub(",", "", sp500$Open))
sp500<-sp500[(nrow(sp500)-2000):nrow(sp500), ]
OpenDiff <- diff(sp500$Open)
plot_me <- OpenDiff


index <- seq(1,length(plot_me), by = 120)
#index <- 1:nrow(sp500)
plot(plot_me, type = "l")

# ARMA MODEL FITTING 
fit <- arima(plot_me[index] , order = c(1, 0,0))
fit
d <- 1

# Get residuals, center them  (centered kernel)
u <- fit$residuals
u <- u - mean(u)
big_T <- length(u)

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
