	ui <- fluidPage(
		radioButtons("dist", "Distribution type:",
					 c("Normal" = "norm",
					   "Uniform" = "unif",
					   "Log-normal" = "lnorm",
					   "Exponential" = "exp"),inline=TRUE),
		plotOutput("distPlot")
	)
	
	server <- function(input, output) {
		output$distPlot <- renderPlot({
			dist <- switch(input$dist,
						   norm = rnorm,
						   unif = runif,
						   lnorm = rlnorm,
						   exp = rexp,
						   rnorm)
			
			hist(dist(500))
		})
	}
	
	shinyApp(ui, server)


	# Not run: 
	library(shiny)
	library(pixiedust)
	library(shinydust)
	
	server <- shinyServer(function(input, output) {
		output$table <-
			renderText({
				cbind(rownames(mtcars), mtcars) %>%
					radioTable(inputId = "chooseCar",
							   label = "",
							   choices = paste0("car", 1:nrow(mtcars)),
							   table_label = "Select a Vehicle",
							   pixie = . %>%
							   	sprinkle(bg_pattern_by = "rows") %>%
							   	sprinkle_table(pad = 7) %>%
							   	sprinkle_colnames("rownames(mtcars)" = "",
							   					  control = ""))
			})
		
		output$choice <- renderText(input$chooseCar)
	})
	
	ui <- shinyUI(fluidPage(
		wellPanel(
			verbatimTextOutput("choice"),
			uiOutput("table")
		)
	))
	
	shinyApp(ui = ui, server = server)
	