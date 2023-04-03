

start_plot <- function(plot_name){
  reso <- 300
  length <- 4
  png(paste(plot_name, ".png", sep = ""),
      units="in", 
      res = reso, 
      width = 6, 
      height = 3)
}
