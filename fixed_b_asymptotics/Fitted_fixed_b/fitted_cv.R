
fitted_cv <- "https://raw.githubusercontent.com/rpkgarcia/DissertationSim/main/fixed_b_asymptotics/Fitted_fixed_b/fitted_CV.csv/"

source(paste(url, "est_autocov.R", sep = ""))
source(paste(url, "kernels.R", sep = ""))


# -------------------------------------------------------------------------
# Fitted CV ---------------------------------------------------------------
# -------------------------------------------------------------------------
fitted_model <- function(d, cv_matrix, alpha){
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
# Get CV ---------------------------------------------------------------
# -------------------------------------------------------------------------

get_cv <- function(d, alpha, the_kernel, is_lugsail, new_b){
  # Read in all fitted values 
  the_fits <- read.csv(fitted_cv)
  
  # Pull out only the values you need 
  index <- the_fits$kernel == the_kernel & the_fits$is_lugsail == is_lugsail &
    the_fits$alpha == alpha & the_fits$d == d
    
  coefficients <- the_fits[index, c("beta1", "beta2","beta3")]
  intercept <-  the_fits[index, c("intercept")]
  
  # Fitted value of b 
  new_b<-data.frame(poly(new_b, 3, raw = T))
  cv_by_b <- apply(new_b, 1, function(x) sum(x*fit$coefficients))
  cv_by_b <- cv_by_b + chisq_cv
  
  return(cv_by_b)
}
