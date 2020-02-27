library(shinydashboard)
source(paste0(wd,"/shiny/files_input.R"))

header <- dashboardHeader(title = "sRNAflow" ,dropdownMenuOutput("messageMenu"))

sidebar <- dashboardSidebar(
    sidebarMenu(
        id = "tabs",
        menuItem("Seq Data Input", tabName="Files", icon = icon("cloud-upload"), selected = TRUE),
        menuItem("Select groups", tabName="Groups", icon = icon("cloud-upload"), selected = TRUE),
        menuItem("Configuration", tabName = "Config", icon = icon("th")),
        # The following menus are just displayed when a sample set has been loaded
        menuItem("Report", tabName = "Report", icon = icon("asterisk"))
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
        menuItem("About", tabName = "About"), textOutput("clock")
    )
)

body <- dashboardBody(
    tabItems(
        # First tab content
        # tabItem(tabName = "dashboard",h2("Work process..."),
        #         fluidRow(
        #             box(plotOutput("distPlot", height = 250)),
        #             
        #             box(
        #                 title = "Controls",
        #                 sliderInput("bins", "Number of observations:", 1, 100, 50)
        #             )
        #         )
        # ),
        tabItem("Files",
                h2("Input files"),
                filesInputUI("datafile")
        ),
        tabItem("Groups",
                h2("Select groups"),
                selectGroupsUI("groups")
        ),
        tabItem(tabName = "Config",
                h2("Configuration")
        ),
        tabItem(tabName = "Report",
                h2("Report")
        ),
        tabItem(
            "About",
            box(width=12,
                HTML(
                    "<h2>sRNAflow</h2>
                                                     
                                                     <p>This tool was developed by ... </p>")),
            br(),
            br(),
            box(width=12,
                title="Session Information",
                collapsible=TRUE,
                collapsed=FALSE,
                verbatimTextOutput("session_info"),
                verbatimTextOutput("session_info1")
            )
        )
    )
)

ui <- dashboardPage(header, sidebar, body)

