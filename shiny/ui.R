library(shinydashboard)
source(file.path(wd,"shiny","files_input.R"))

header <- dashboardHeader(title = "sRNAflow" ,dropdownMenuOutput("messageMenu"))

sidebar <- dashboardSidebar(
    sidebarMenu(
        id = "tabs",
        menuItem("Seq Data Input", tabName="Files", icon = icon("cloud-upload"), selected = TRUE),
        menuItem("Select groups", tabName="Groups", icon = icon("tasks")),
        menuItem("Analysis", tabName = "Config", icon = icon("sliders-h")),
        # The following menus are just displayed when a sample set has been loaded
        menuItem("Reports", tabName = "Reports", icon = icon("chart-bar")),
        menuItem("Setup", tabName = "Setup", icon = icon("tools"))
    ),
    div(class="hide_when_sidebar_collapsed",
        # tags$p(class="sidebartext", style="padding-left: 10px;color: #b8c7ce; ", "To start  data, upload a dataset in the 'Seq Data Input' tab."),
        uiOutput("bookmarkBtnUI")
    ),
    # conditionalPanel(
    #     condition = "bins == 3",
    #     sidebarMenu(
    #         menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard"))
    #     )
    # ),
    sidebarMenu(
        menuItem("About", tabName = "About", icon = icon("info")), textOutput("clock")
    )
)

body <- dashboardBody(
    tabItems(
        tabItem("Files",
                h2("Input files"),
                filesInputUI("datafile")
        ),
        tabItem("Groups",
                h2("Select groups"),
                selectGroupsUI("groups"),
                tableOutput("sel")
        ),
        tabItem(tabName = "Config",
                h2("Configuration"),
                fluidRow(
                    column(3,textInput("Exp",   "Experiment ID:", value = Exp, width ='100%')),
                    column(3,selectInput('strategy', 'Strategy:', c("successively","metagenome"), selected = strategy, width ='100%')),
                    column(3,textInput("email",   "Send results to (email):", value = email, width ='100%')),
                    column(3,textInput("smtpServer",   "SMTP Server:", value = smtpServer, width ='100%'))
                ),
                hr(),
                fluidRow(
                    column(4,
                        h3("Trimming"),
                        sliderInput('qc', 'QC trimming:', min = 1,  max = 30,  value = qc, width ='100%')
                    ),
                    column(4,
                        h3("BLAST"),
                        sliderInput('tsize', 'Representative subset size (use ~200 for remote db) :', min = 200,  max = 2000,  value = tsize, width ='100%')
                    ),
                    column(4,autoWidth = TRUE,
                        h3("Differential expression"),
                        fluidRow(
                            column(6,sliderInput('lim', 'Keep hits that more then:', min = 0,  max = 100,  value = lim, width ='100%')),
                            column(6,sliderInput('limS','in more then "x" samples; x=:', min = 0,  max = 100,  value = limS, width ='100%'))
                        )
                    )
                ),
                hr(),
                fluidRow(
                    column(4,selectizeInput('ad3', "3' adapter:", choices = c('TGGAATTCTCGGGTGCCAAGG #Illumina TruSeq Small RNA','ATCACCGACTGCCCATAGAGAG # Ion Torrent','AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCGGCTTG # Universal Illumina adapter',' # Not remove'),options = list(create = TRUE), selected = ad3, width ='100%')),
                    column(4,autoWidth = TRUE,sliderInput('Rep', 'Subsets number:', min = 1,  max = 10,  value = Rep, width ='100%')),
                    column(4,numericInput('log2FoldChange', 'log2FoldChange theshold:', log2FoldChange, min = 0, max = 100, width ='100%'))
                ),
                hr(),
                fluidRow(
                    column(4,selectizeInput('ad5', "5' adapter:", choices = c('GTTCAGAGTTCTACAGTCCGACGATC # Illumina TruSeq Small RNA','CCAAGGCG # Ion Torrent',' # Not remove'),options = list(create = TRUE), selected = ad5, width ='100%')),
                    column(4,selectInput('specie', 'Main specie:', species, selected = specie, width ='100%')),
                    column(4,numericInput('padj', 'Adjusted p-value theshold:', padj, min = 0, max = 1, width ='100%'))
                ),
                hr(),
                fluidRow(
                    column(4,sliderInput("sizerange", "Size range", min = 10,max = 300, value = sizerange, width ='100%')),
                    column(4,selectInput('blast', 'BLAST to (nr/nt will be used if absent local db):', c("main specie & bacteria+","nr/nt"), selected = blast, width ='100%')),
                    column(4,br(),br(),actionButton("start", "Start analysis",icon = icon("rocket"), width ='100%'))
                ),
                verbatimTextOutput("Config")
        ),
        tabItem(tabName = "Reports",
                h2("Reports"),
                hr(),
                showReports("reports"),
                actionButton("refresh_reports", "Refresh reports list",icon = icon("sync"), width ='100%')
        ),
        tabItem(tabName = "Setup",
                h2("Setup"),
                hr(),
                actionButton("blastdb", "Create/Update local BLAST db",icon = icon("sync"), width ='100%'),
                verbatimTextOutput("blastdb"),
                br(),hr(),
                actionButton("gtfdb", "Update GTF files",icon = icon("sync"), width ='100%'),
                verbatimTextOutput("gtfdb"),
                br(),hr(),
                verbatimTextOutput("Setup")
        ),
        tabItem(
            "About",
            box(width=12,
                HTML(
                    "<h2>sRNAflow</h2>
                    <a href=https://github.com/zajakin/sRNAflow>https://github.com/zajakin/sRNAflow</a>
                    <p>This tool was developed by Pawel Zayakin.</p>")),
            br(),
            box(width=12,
                title="Session Information",
                collapsible=TRUE,
                collapsed=TRUE,
                verbatimTextOutput("session_info")
            )
        )
    )
)

ui <- dashboardPage(header, sidebar, body)

