# install.packages("remotes")
# remotes::install_github("fbreitwieser/shinyFileTree")
library(shinyFileTree)

uploadFilePanel <- function(ns) {
	tabPanel("Upload files",
		DT::dataTableOutput("filesUploaded"),
		fileInput(
		 	"file_upload",
		 	label="",
		 	placehold = "Upload files",
		 	width = "100%",
		 	multiple = TRUE
		)
	)
}

serverDataPanel <- function(ns) {
	tabPanel("Use data on server",
		{if (!dir.exists(paste0(wd,"/data/input"))) dir.create(paste0(wd,"/data/input"))
			DT::dataTableOutput('serverFiles')
		} ,
		shinyFileTree::shinyFileTreeOutput(ns('file_tree'))
	)
}

exampleDataPanel <- function(ns) {
	tabPanel("Example samples",
		{if (!dir.exists(paste0(wd,"/data/example-samples"))) dir.create(paste0(wd,"/data/example-samples"))
			DT::dataTableOutput('examples')
		}
	)
}

filesUploadRow <- function(x){
	return(c(file=gsub(" ","_",paste0("data/input/",x$name)),
			size=humanReadable(x$size),
	  		date=format(Sys.time(),"%d.%m.%Y %H:%M:%OS")))
}

filesInputUI <- function(id) {
	ns <- NS(id)
	shiny::tagList({
		shinydashboard::tabBox(
			width = 12,
			title = "Data Source",
			uploadFilePanel(ns),
			serverDataPanel(ns),
			exampleDataPanel(ns)
		)
	},
	hr(),
	h2("Selected files:"),
	DT::dataTableOutput("filesIn"),
	# fluidRow(
	# 	column(9, DT::dataTableOutput("filesIn")),
	# 	column(3, verbatimTextOutput("filesIn_selected"))
	# ),

	# h1('A Client-side Table'),
	# fluidRow(
	# 	column(6, DT::dataTableOutput('x1')),
	# 	column(6, plotOutput('x2', height = 500))
	# ),
	# 
	# h1('A Server-side Table'),
	# fluidRow(
	# 	column(9, DT::dataTableOutput('x3')),
	# 	column(3, verbatimTextOutput('x4'))
	# )
	)
}
selectGroupsUI <- function(id) {
	ns <- NS(id)
	shiny::tagList(
		hr(),
		h2("Selected files:")
		# DT::dataTableOutput("filesIn")
	)
}
