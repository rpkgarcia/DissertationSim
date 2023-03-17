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


