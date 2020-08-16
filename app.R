#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
options(shiny.maxRequestSize = 30000*1024^2)
source("shiny/global.R")
source("shiny/ui.R")
source("shiny/server.R")
#TODO Upload big files
server_options<-list(port=4444,shiny.maxRequestSize=10000*1024^2)
# timeout<-2147483
# server_options<-append(server_options,c(app_init_timeout=timeout,shiny.app_init_timeout=timeout,
# shiny.app_idle_timeout=timeout,app_connection_timeout=timeout,shiny.app_connection_timeout=timeout,
# app_read_timeout=timeout,shiny.app_read_timeout=timeout,app_idle_timeout=timeout,shiny.trace=T))

# Run the application 
shinyApp(ui = ui, server = server, options = server_options)
# setwd("shiny"); shiny::runApp(display.mode="showcase")
