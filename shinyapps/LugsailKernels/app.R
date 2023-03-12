# This app is to illustrate the difference between the different lugsails

library(stringr)

# Support Functions -------------------------------------------------------
url <- "https://raw.githubusercontent.com/rpkgarcia/DissertationSim/main/functions/"

source(paste(url, "kernels.R", sep = ""))


# Shiny App ---------------------------------------------------------------


ui <- fluidPage(
  
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      sliderInput("big_T_value", "Select Sample Size:", 
                  value = 200,
                  min=100, 
                  max = 100000, 
                  step = 100),
      sliderInput("b_value", "Select Bandwidth:", 
                  value = .071,
                  min = .02,
                  max = 1, 
                  step = .02),
      selectInput(inputId = "kernel", 
                  label = "Select Kernel Function:", 
                  choices = c("Bartlett" = "bartlett", 
                              "Parzen" = "parzen", 
                              "QS" = "qs")),  
      selectInput("x_range", 
                  "X-Axis Range", 
                  c("Full View", "Half View")),
      checkboxGroupInput("kernelOptions", 
                         "Select Kernels",
                         c("Zero", "Adaptive", "Over"))
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot")
    )
  )
  
)

server <- function(input, output, session) {
  # observe({
  #   x <- input$lugsailOptions
  #   
  #   # Can use character(0) to remove all choices
  #   if (is.null(x))
  #     x <- character(0)
  #   
  #   # Can also set the label and select items
  #   updateSelectInput(session, "inSelect",
  #                     label = paste("Select input label", length(x)),
  #                     choices = x,
  #                     selected = tail(x, 1)
  #   )
  # })
  
  
  output$distPlot <- renderPlot({
    
    # generate bins based on input$bins from ui.R
    if(input$x_range == "Full View"){
      x = seq(-1.2, 1.2, length.out = 2000)
    } else{
      x = seq(0, 1.2, length.out = 2000)
    }
    
    y = sapply(x, input$kernel)
    
    if(input$kernel == "bartlett"){
      mother <- "Bartlett"
      q = 1
    } else { 
      q = 2
      if(input$kernel == "parzen"){
        mother <- "Parzen"
      } else {
        mother <- "QS"
      }
    }
    
    # Get Maximum values for y-axis
    y_max = c(1)
    if(any(str_detect(input$kernelOptions, "Adaptive"))){
      lug_para <- get_lugsail_parameters(big_T = input$big_T_value, 
                                         b = input$b_value,
                                         q = q,
                                         method ="Adaptive")
      y_lugsail = sapply(x,lugsail,
                         lugsail_parameters = lug_para,
                         the_kernel = eval(parse(text=input$kernel)))
      y_max = c(y_max, max(y_lugsail))
    } 
    
    if(any(str_detect(input$kernelOptions, "Over"))){
      lug_para <- get_lugsail_parameters(method ="Over", q = q)
      y_lugsail = sapply(x,lugsail,
                         lugsail_parameters = lug_para, 
                         the_kernel = eval(parse(text=input$kernel)))
      y_max = c(y_max, max(y_lugsail))
    } 
    
    plot(x, y, xlab = "x",
         ylab = "weight", 
         type = "l", ylim = c(0, max(y_max)),
         lwd = 2)
    
    # Plot extra curves as needed
    the_colors = rainbow(3)
    
    keep = which(c("Zero","Adaptive", "Over") %in% input$kernelOptions)
    the_colors = the_colors[keep]
    the_lwd = 2:4
    the_lwd = the_lwd[keep]
    
    if(length(input$kernelOptions) != 0 ){
      for(i in 1:length(input$kernelOptions)){
        
        option = input$kernelOptions[i]
        lug_para <- get_lugsail_parameters(big_T = input$big_T_value, 
                                           b = input$b_value, 
                                           q = q,
                                           method = option)
        y_lugsail <- sapply(x,lugsail,
                            lugsail = lug_para,
                            the_kernel = eval(parse(text=input$kernel)))
        
        lines(x, y_lugsail, col = the_colors[i], lwd = 2, lty = the_lwd[i])
      }
      
    }
    
    
    legend("topright", 
           c(mother, input$kernelOptions), 
           lty = c(1, the_lwd), 
           lwd = 2, 
           col = c("black", the_colors))
  })
  
}

shinyApp(ui, server)
