setwd("~/Documents/GitHub/DissertationSim/fixed_b_asymptotics/1_dimensional/psd_issues")
try_b <-  seq(0.005, .99, by = 0.005)
the_y_lim <- c(0, .3)
the_x_lim <- c(0, .50)

psd_bartlett <- read.csv("psd_issues/bartlett.csv", row.names=1)
psd_parzen <- read.csv("psd_issues/parzen.csv", row.names=1)
psd_qs <- read.csv("psd_issues/qs.csv", row.names=1)


plot_psd_rate <- function(psd_counts, the_ylab = "", the_main){
  psd_zero <- psd_counts[[1]]
  psd_adapt <- psd_counts[[2]]
  psd_over <- psd_counts[[3]]
  #the_main <- gsub("psd_", "", deparse(substitute(psd_counts)))
  #the_main <- StrCap(the_main, method = "first")
  
  plot(try_b, psd_zero, 
       xlab = "b", 
       ylab = the_ylab, 
       col = "red", type = "l", lty = 2,  
       xlim = the_x_lim, 
       ylim = the_y_lim, 
       main = the_main)
  lines(try_b, psd_adapt, col = "blue", lty = 3)
  lines(try_b, psd_over, col = "green", lty = 4)
}



start_plot <- function(plot_name){
  reso <- 300
  length <- 4
  png(paste(plot_name, ".png", sep = ""),
      units="in", 
      res = reso, 
      width = 6, 
      height = 3)
}

start_plot("psd_proportion")
par(mfrow = c(1, 3))
plot_psd_rate(psd_bartlett, "Proportion Corrected", "Bartlett")
plot_psd_rate(psd_parzen, the_main="Parzen")
plot_psd_rate(psd_qs, the_main="Quadratic Spectral")
dev.off()
