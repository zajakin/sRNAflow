#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

source("shiny/global.R")
source("shiny/ui.R")
source("shiny/server.R")

# Run the application 
shinyApp(ui = ui, server = server, options = list(port=4444))
# setwd("shiny"); shiny::runApp(display.mode="showcase")
