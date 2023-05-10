
try_b <-  seq(0.005, .99, by = 0.005)
the_y_lim <- c(0, .03)
the_x_lim <- c(0, .15)

psd_bartlett <- read.csv("psd_issues/bartlett.csv", row.names=1)
psd_parzen <- read.csv("psd_issues/parzen.csv", row.names=1)
psd_qs <- read.csv("psd_issues/qs.csv", row.names=1)


plot_psd_rate <- function(psd_counts, the_ylab = ""){
  psd_zero <- psd_counts[[1]]
  psd_adapt <- psd_counts[[2]]
  psd_over <- psd_counts[[3]]
  the_main <- gsub("psd_", "", deparse(substitute(psd_counts)))
  #the_main <- StrCap(the_main, method = "first")
  
  plot(try_b, psd_zero, 
       xlab = "b", 
       ylab = the_ylab, 
       col = "blue", type = "l", 
       xlim = the_x_lim, 
       ylim = the_y_lim, 
       main = the_main)
  lines(try_b, psd_adapt, col = "purple")
  lines(try_b, psd_over, col = "red")
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
plot_psd_rate(psd_bartlett, "Proportion Corrected")
plot_psd_rate(psd_parzen)
plot_psd_rate(psd_qs)
dev.off()
