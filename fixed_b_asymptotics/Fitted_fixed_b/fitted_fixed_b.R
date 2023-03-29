# Create a table with the fixed-b coefficients 
#         CV_alpha(b) = beta_0 + b*beta_1 + b^2*beta_2

# Coefficents are specific to b, and kernel 



# -------------------------------------------------------------------------
# Load Functions ----------------------------------------------------------
# -------------------------------------------------------------------------
library(stringr)
library(dplyr)


main <- "https://raw.githubusercontent.com/rpkgarcia/DissertationSim/main/"
url <- paste(main, "fixed_b_asymptotics/Fitted_fixed_b", sep = "") 
source(paste(url, "fitted_cv.R", sep = ""))


# -------------------------------------------------------------------------
# Load Functions ----------------------------------------------------------
# -------------------------------------------------------------------------

# Get a list of all the CV files 
setwd("~/Documents/GitHub/DissertationSim/fixed_b_asymptotics/Fixed_b_CV_Tables")
files <- list.files()

# Read in the file, and create a fitted linear regression line 
# Keep information about the file: mother kernel, lugsail, alpha, d

the_fits <- lapply(files, function(the_file){
  cv_matrix <- read.csv(the_file)
  the_kernel <- gsub("_","", str_extract(the_file, "[[:alpha:]]+_"))
  is_lugsail <- grepl("Lugsail", the_file)
  alpha <- as.numeric(paste(".", gsub("\\.csv","",str_extract(the_file, "[[:digit:]]+\\.csv")), sep=""))

  
  
  the_fits <- sapply(1:(ncol(cv_matrix)-1), function(d){
    fit <- fitted_model(d, cv_matrix, alpha)
    fit <- c(fit, kernel = the_kernel ,is_lugsail = is_lugsail, alpha = alpha)
  })
  the_fits <- data.frame(t(the_fits))
})

# Store the results for easy access 
setwd("../Fitted_fixed_b")
the_fits <- bind_rows(the_fits, .id = "file")
write.csv(the_fits, "fitted_CV.csv",row.names = F) 

