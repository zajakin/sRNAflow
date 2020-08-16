#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
timeout<-2147483
options(shiny.maxRequestSize = 30000*1024^2,port=4444,app_init_timeout=timeout,shiny.app_init_timeout=timeout,
        shiny.app_idle_timeout=timeout,app_connection_timeout=timeout,shiny.app_connection_timeout=timeout,
        app_read_timeout=timeout,shiny.app_read_timeout=timeout,app_idle_timeout=timeout,shiny.trace=T)
source("shiny/global.R")
source("shiny/ui.R")
source("shiny/server.R")

# Run the application 
shinyApp(ui = ui, server = server)
