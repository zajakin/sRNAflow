# devtools::install_github("fbreitwieser/shinyFileTree", type = "source")
# library(shinyFileTree)

uploadFilePanel <- function(ns) {
	tabPanel("Upload files",
		DT::dataTableOutput("filesUploaded"),
		fileInput(
		 	"file_upload",
		 	label="",
		 	placehold = "Upload files",
		 	width = "100%",
		 	multiple = TRUE
		),
		uiOutput("filesUpload")
	)
}

serverDataPanel <- function(ns) {
	tabPanel("Use data on server",
		{	if(!dir.exists(file.path(wd,"www","upload"))) dir.create(file.path(wd,"www","upload"),recursive = TRUE, mode = "0777")
			DT::dataTableOutput('serverFiles')
		},
		actionButton("refresh_serverFiles", "Refresh files list",icon = icon("sync"), width ='100%')
		# , shinyFileTree::shinyFileTreeOutput(ns('file_tree'))
	)
}

exampleDataPanel <- function(ns) {
	tabPanel("Example samples",
		{ if(!dir.exists(file.path(wd,"www","upload","example-samples"))) dir.create(file.path(wd,"www","upload","example-samples"),recursive = TRUE, mode = "0777")
			DT::dataTableOutput('examples')
		},
		fluidRow(
			column(4,actionButton("download_examples", "Download 39 EV examples",icon = icon("vials"), width ='100%')),
			column(4,actionButton("download_examples2", "Download 40 cells examples",icon = icon("vials"), width ='100%')),
			column(4,actionButton("refresh_examples", "Refresh examples list",icon = icon("sync"), width ='100%'))
		)
	)
}

filesUploadRow <- function(x){
	return(c(file=gsub(" ","_",x$name),
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
	    actionButton("clear_filesIn", "Clear list of selected files",icon = icon("trash"), width ='100%')
		# ,fluidRow(
		# 	column(9, DT::dataTableOutput("filesIn")),
		# 	column(3, verbatimTextOutput("filesIn_selected")),
		# 	column(6, plotOutput('x2', height = 500))
		# )
	)
}

selectGroupsUI <- function(id) {
	ns <- NS(id)
	shiny::tagList(
		hr(),
		DT::dataTableOutput("groups")
	)
}

showReports<- function(id) {
	ns <- NS(id)
	shiny::tagList(
		hr(),
		DT::dataTableOutput("reports")
	)
}
