# -------------------------------------------------------------------------
# Load Functions ----------------------------------------------------------
# -------------------------------------------------------------------------
url <- "https://raw.githubusercontent.com/rpkgarcia/DissertationSim/main/simulations/"
source(paste(url, "multi_sim_f_stat.R", sep = ""))
source(paste(url, "fitted_cv.R", sep = ""))


# -------------------------------------------------------------------------
# Simulation Parameters (Adjust me) ---------------------------------------
# -------------------------------------------------------------------------

# How many replicates
# KV005 used 50,0000
num_replicates <- 100

# Sample size of each replicate
big_T <- 200

# Maximum number of dimensions
d_max <- 2

# Correlation strength 
rho_y <- 0.7

# optimal b: Bartlett
opt_b <- mapply(optimal_b, rho_y, alpha = 0.05, big_T=big_T, kernel=1, 
                mother=2, d=1:d_max) 

#  What bs to use 
try_b <-  c(0.010, 0.020, 0.025, 
            0.030, 0.035, 0.040, 0.045, 0.050, 0.055, 
            0.065, 0.075, 0.080, 0.092, 0.095, 
            0.100, 0.120, 0.140, 0.160, 0.180, 0.200, 
            0.220, 0.240, 0.260, 0.280, 0.300, 0.350,
            0.400, 0.450, 0.500, 0.555, 0.600, 0.650, 
            0.700, 0.750, 0.995, opt_b)
# Small for now
try_b <- c(0.010, 0.030,  0.050, 0.07, 
           0.10, 0.20, 0.30, 0.40, 0.5, opt_b)
# -------------------------------------------------------------------------
# Run Simulation  ---------------------------------------------------------
# -------------------------------------------------------------------------
all_F_stats = replicate(num_replicates, 
                        simulate_f_stat(big_T = big_T, d_max = d_max))

# [ #,  ,  ,  ] : rows, the different b values [1:length(try_b)]
# [  , #,  ,  ] : columns, the different d values [1:max_d]
# [  ,  , #,  ] : the original kernel (1), and lugsail (2) [1:2]
# [  ,  ,  , #] : each simulation [1:num_replicates]

dimnames(all_F_stats)[[1]] = try_b
dimnames(all_F_stats)[[2]] = 1:d_max 
dimnames(all_F_stats)[[3]] = c("Bartlett", "Bartlett_Lugsail") 
dimnames(all_F_stats)[[4]] = paste("sim", 1:num_replicates, sep="")

# -------------------------------------------------------------------------
# Get CV 
# -------------------------------------------------------------------------

alpha <- 0.05 
kernel <- "Bartlett"

file <- "../fixed_b_asymptotics/Fixed_b_CV_Tables/Bartlett_05.csv"
the_cv <- read.csv(file)

CV <- sapply(1:d_max, fitted_model, cv_matrix=the_cv, new_b=try_b, alpha=0.05)
colnames(CV) <- 1:d_max
rownames(CV) <- try_b



# -------------------------------------------------------------------------
# Size --------------------------------------------------------------------
# -------------------------------------------------------------------------


size <- array(0, dim = dim(all_F_stats)[1:3], 
              dimnames = list(try_b, 
                              1:d_max, 
                              c("Bartlett", "Bartlett_Lugsail")))


for(index_b in 1:length(try_b)){
  for(index_d in 1:d_max){
    total <- sum(all_F_stats[index_b, index_d, "Bartlett_Lugsail", ] > CV[index_b, index_d])
    size[index_b, index_d, "Bartlett_Lugsail"] <- total/num_replicates
  }
}

size_df <- data.frame(b = rep(try_b, times = d_max), 
                      d = rep(1:d_max, each = length(try_b)), 
                      size = c(size))
setwd("saved_stats")
write.csv(size_df, "size.csv")
setwd("..")

plot(try_b, size[, 1, 2], ylim = c(0, .10))
abline(h = 0.05, lty =2, col = "Red")
opt_b_index <- (length(try_b)-1)
points(try_b[opt_b_index], size[opt_b_index , 1, 2], pch = 8, col = "green")


# -------------------------------------------------------------------------
# Power  --------------------------------------------------------------------
# -------------------------------------------------------------------------
ncp <- seq(0, 2.5, by = .25)
ncp <- c(ncp, 3:5)

power <- array(0, dim = c(dim(all_F_stats)[1:3], length(ncp)), 
               dimnames = list(try_b, 
                               1:d_max, 
                               c("Bartlett", "BartlettLugsail"), 
                               ncp))

the_seed <- 2

for(index_ncp in 1:length(ncp)){
  set.seed(the_seed)
  mean_vec <- rep(ncp[index_ncp], d_max)
  all_F_stats = replicate(num_replicates, 
                          simulate_f_stat(big_T, d_max, mean_vec))
  for(index_b in 1:length(try_b)){
    for(index_d in 1:d_max){
      print(paste("ncp=", index_ncp, "; try_b=", index_b, "; index_d=", index_d, sep = ""))
      total <- sum(all_F_stats[index_b, index_d, 2, ] < CV[index_b, index_d])
      power[index_b, index_d, 2, index_ncp] <- total/num_replicates
    }
  }
  
}

power_df <- data.frame(b = rep(try_b, times = d_max*length(ncp)), 
                       d = rep(1:d_max, each = length(try_b), times = length(ncp)),
                       ncp = rep(ncp, each = length(try_b)*d_max), 
                       power = c(power))
setwd("saved_stats")
write.csv(size_df, "power.csv")
setwd("..")

opt_b_index <- (length(try_b)-1)
plot(ncp, power[opt_b_index, 1, 2, ], ylim = c(0, 1), type = "b")
abline(h = 0.05, lty = 3)
