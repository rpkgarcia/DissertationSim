# This is only for one-dimension 
library("stringr")
library(dplyr)
library(tidyverse)


# -------------------------------------------------------------------------
# Fitted CV ---------------------------------------------------------------
# -------------------------------------------------------------------------
fitted_model <- function(d=1, cv_matrix, alpha){
  try_b <- cv_matrix[,1]
  alpha_levels <- 1-alpha 
  chisq_cv <-  qchisq(alpha_levels, df = d)/d 
  specific_cvs <- cv_matrix[,d+1] - chisq_cv
  
  # Get a CV fitted as a function of b 
  fit <-  lm(specific_cvs ~ 0+  poly(try_b, 3, raw = T)) 
  fit <- summary(fit)
  
  return_me <- c(chisq_cv, fit$coefficients[,1], fit$adj.r.squared, d)
  names(return_me) <- c("intercept", "beta1", "beta2", "beta3", "adj.r.sq", "d")
  
  return(return_me)
}


# -------------------------------------------------------------------------
# Fitted CV ---------------------------------------------------------------
# -------------------------------------------------------------------------

files <- list.files()
files <- grep("\\.csv",files, value = T)
files <- files[-which(files == "fitted_CV.csv")]

# Read in the file, and create a fitted linear regression line 
# Keep information about the file: mother kernel, lugsail, alpha, d

the_fits <- data.frame(NA, NA, NA,NA, NA, NA, NA, NA, NA)
colnames(the_fits) <- c("intercept", "beta1", "beta2", "beta3", "adj.r.sq","d",
                        "kernel", "lugsail", "alpha")
for(the_file in files){
  cv_matrix <- read.csv(the_file)
  the_file <- strsplit(the_file, "_")[[1]]
  mother_kernel <- the_file[1]
  lugsail_type  <- the_file[2]
  alpha <- as.numeric(paste(".", gsub("\\.csv","",the_file[3]), sep=""))
  fit <- fitted_model(1, cv_matrix, alpha)
  fit <- c(fit, kernel = mother_kernel , lugsail = lugsail_type, alpha = alpha)
  the_fits <- rbind(the_fits, fit)
}

# Store the results for easy access 
the_fits <- the_fits[-1, ]
write.csv(the_fits, "fitted_CV.csv",row.names = F) 
