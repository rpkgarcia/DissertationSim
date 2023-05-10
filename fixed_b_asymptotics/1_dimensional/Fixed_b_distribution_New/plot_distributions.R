# Plot Distribution 

index_b <- 15 # 0.07 
the_x_lim <- c(0,6)
the_y_lim = c(0, .9)

index <- 1:5000

start_plot <- function(plot_name){
  reso <- 300
  length <- 4
  png(paste(plot_name, ".png", sep = ""),
      units="in", 
      res = reso, 
      width = 6, 
      height = 4)
}


plot_density <- function(the_kernel, the_ylab = ""){
  start_plot(the_kernel)

  Mother <- read.csv(paste("Fixed_b_distribution_New/", the_kernel, "_Mother",
                           ".csv", sep = ""))
  Zero <- read.csv(paste("Fixed_b_distribution_New/", the_kernel, "_Zero",
                         ".csv", sep = ""))
  Adapt <- read.csv(paste("Fixed_b_distribution_New/", the_kernel, "_Over",
                          ".csv", sep = ""))
  Over <- read.csv(paste("Fixed_b_distribution_New/", the_kernel, "_Adapt",
                         ".csv", sep = ""))

  
  plot(density(Mother[index,index_b], from = 0, to = 10), xlim = the_x_lim, ylim = the_y_lim,
       main = paste(the_kernel, "~new~",max(index)), ylab = the_ylab, xlab  = "F-Statistic")
  curve(dchisq(x, 1), add = T, lty = 2, col = "grey", lwd = 2)
  lines(density(Zero[index,index_b], from = 0, to = 10), xlim = the_x_lim, col = "red", lty = 2)
  lines(density(Adapt[index,index_b], from = 0, to = 10), xlim = the_x_lim, col =  "blue", lty = 3)
  lines(density(Over[index,index_b], from = 0, to = 10), xlim = the_x_lim, col = "green", lty =4)
  dev.off()
  
}

#start_plot("density")
# par(mfrow= c(1, 1))
plot_density("Bartlett", "Density")
plot_density("Parzen", "Density")
plot_density("QS", "Density")
#dev.off()


# library(ggplot2)
# plot_me <- rbind(data.frame(F_stat =Mother[,index_b], Kernel = "Mother"), 
#                  data.frame(F_stat =Zero[,index_b], Kernel = "Zero"),
#                  data.frame(F_stat =Adapt[,index_b], Kernel = "Adapt"),
#                  data.frame(F_stat =Over[,index_b], Kernel = "Over")
#                  )
# 
# the_seq <- seq(0, 10, length.out = 200)
# ggplot(plot_me, aes(F_stat))+
#   geom_density(aes(color = Kernel, linetype = Kernel))+
#   stat_function(fun = dchisq, args = list(df = 1), 
#                 color = "grey", lty = 2)+ 
#   coord_cartesian(xlim= c(0, 10),ylim = c(0, 2))+ 
#   theme_classic()+ 
#   scale_color_manual(values = c("Black", "Red",  "purple", "Blue", "Grey"), 
#                      limits = c("Mother", "Zero", "Adapt","Over", "Small-b"))+
#   scale_linetype_manual(values = c(1, 2,  3, 4, 2), 
#                      limits = c("Mother", "Zero", "Adapt","Over", "Small-b"))+
#   labs(x = "F-Statistics", y = "Density")
#   
