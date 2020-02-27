# autoInvalidate <- reactiveTimer(10000)
# observe({
#     autoInvalidate()
#     cat(".")
# })

server <- function(input, output) {
    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        x    <- faithful[, 2]
        bins <- seq(min(x), max(x), length.out = input$bins + 1)
        
        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'white')
    })
    output$messageMenu <- renderMenu({
        # Code to generate each of the messageItems here, in a list. This assumes
        # that messageData is a data frame with two columns, 'from' and 'message'.
        messageData<-rbind(c(from="",message=""))
        msgs <- apply(messageData, 1, function(row) {
            messageItem(from = row[["from"]], message = row[["message"]])
        })
        
        dropdownMenu(type = "messages", .list = msgs)
        dropdownMenu(type = "notifications", .list = msgs)
        dropdownMenu(type = "tasks", .list = msgs)
    })
    output$clock <- renderText({
        invalidateLater(1000)
        format(Sys.time(),"%d.%m.%Y %H:%M:%OS")
    })
    output$filesUpload <- renderText({
        req(input$file_upload)
        if(!dir.exists(paste0(wd,"/data/input"))) dir.create(paste0(wd,"/data/input"),recursive = TRUE)
        out<-gsub(" ","_",paste0(wd,"/data/input/",input$file_upload$name))
        file.copy(input$file_upload$datapath,out)
        file.remove(input$file_upload$datapath)
        out
    })
    # output$contents <- renderTable({
    #     inFile <- input$file_upload
    #     if (is.null(inFile))
    #         return(NULL)
    #     return(input$file_upload)
    # })
    output$filesUploaded = DT::renderDataTable({
        # req(input$file_upload)
        # if (length(input$filesIn_rows_selected)) observeEvent(input$filesUploaded_rows_selected,shinyjs::reset("filesUploaded_rows_selected"))
        # if ((!is.null(input$file_upload) | length(input$file_upload)>0) & length(input$filesIn_rows_selected)==0){ 
        if (!is.null(input$file_upload)) filesUploaded<<-rbind(filesUploaded,filesUploadRow(input$file_upload))
        # observeEvent(input$file_upload,shinyjs::reset("file_upload"))
        filesUploaded
    }, server = TRUE)

    output$examples = DT::renderDataTable({
        if (length(input$filesIn_rows_selected)) observeEvent(input$examples_rows_selected,shinyjs::reset("examples_rows_selected"))
        examples
    }, server = TRUE)
    
    output$serverFiles = DT::renderDataTable({
        if (length(input$filesIn_rows_selected)) observeEvent(input$serverFiles_rows_selected,shinyjs::reset("serverFiles_rows_selected"))
        serverFiles
    }, server = TRUE)
    
aaa<<-input
    output$filesIn = DT::renderDataTable({
    # output$filesIn = renderPrint({
        # cat(input)
        # req(input$examples_row_selected)
        # examples<-dir(path = paste0(wd,"example-samples"),full.names = FALSE, recursive = TRUE, include.dirs = TRUE)
        # examples<-cbind(examples,humanReadable(file.info(paste0(wd,"/example-samples/",examples))$size),format(file.info(paste0(wd,"/example-samples/",examples))$mtime,"%d.%m.%Y %H:%M:%OS"))
        # rbind()
        # examples
        s = input$filesIn_rows_selected
        if (length(s)) filesIn <- rbind(filesIn[-s,])
        s = input$filesUploaded_rows_selected
        # input$filesUploaded_rows_selected <- NULL
        observeEvent(input$filesUploaded_rows_selected,shinyjs::reset("filesUploaded_rows_selected"))
        if (length(s)) filesIn <- rbind(filesUploaded[s,],filesIn)
        s = input$examples_rows_selected
        # input$examples_rows_selected <- NULL
        observeEvent(input$examples_rows_selected,shinyjs::reset("examples_rows_selected"))
        if (length(s)) filesIn <- rbind(examples[s,],filesIn)
        s = input$serverFiles_rows_selected
        # input$serverFiles_rows_selected <- NULL
        observeEvent(input$serverFiles_rows_selected,shinyjs::reset("serverFiles_rows_selected"))
        if (length(s)) filesIn <- rbind(serverFiles[s,],filesIn)
        # colnames(filesIn)<-colnames(examples)
        if(nrow(filesIn)>0) colnames(filesIn)<-c("file","size","date")
        filesIn <- unique(filesIn)
        filesIn <<- filesIn
        # filesIn$dependencies <-
        filesIn
    }, server = TRUE)

# output$groups = DT::renderDataTable({
output$groups = DT::renderDataTable({
        if(length(c(input$filesUploaded_rows_selected,input$examples_rows_selected,input$serverFiles_rows_selected,input$filesIn_rows_selected)))
            groups <- cbind(rbind(filesIn),test="",control="",ignore="")
        colnames(groups)<-c("file","size","date","test","control","ignore")
        # groups<-list()
        if(nrow(groups)>0){
            rownames(groups)<-paste0("S",1:nrow(groups))
            groups[,"test"]   <-paste0('<input type="radio" name="',rownames(groups),'" value="1" //>')
            groups[,"control"]<-paste0('<input type="radio" name="',rownames(groups),'" value="0" //>')
            groups[,"ignore"] <-paste0('<input type="radio" name="',rownames(groups),'" value="-1" //>')
            # for(s in 1:nrow(groups)){
            #     groups[s,"test"]<-radioButtons(paste0("S",s),"Distribution type:",
            #                     choiceNames = list(HTML("<td>Norm</td>"),HTML("<td>Control</td>"),
            #                                        HTML("<td>Ignore</td>")),
            #                                    choiceValues = list("test", "control", "ignore"),inline=TRUE)
            # }
        } #else filesIn<-matrix(0, ncol=5, nrow=0,dimnames = list(month.abb, LETTERS[1:5]))
        groups <<- groups
        # print(groups)
        groups
    }, server = TRUE, escape = FALSE, selection = 'none', options = list(dom = 't', paging = FALSE, ordering = FALSE), 
    callback = DT::JS("table.rows().every(function(i, tab, row) {
          var $this = $(this.node());
          $this.attr('id', this.data()[0]);
          $this.addClass('shiny-input-radiogroup');
        });
        Shiny.unbindAll(table.table().node());
        Shiny.bindAll(table.table().node());")
)  #             sel[1:nrow(filesIn)]

output$sel = renderPrint({
    print(length(input))
    if(length(c(input$filesUploaded_rows_selected,input$examples_rows_selected,input$serverFiles_rows_selected,input$filesIn_rows_selected)))
        print(rownames(groups))
    # if(nrow(groups)>0 & length(input$S1)>0)
        # if(length(mget(paste0("input$",rownames(groups))))>0){
            # print(rownames(groups))
            # sel<-str(sapply(rownames(groups), function(i) input[[i]]))
            # sel<-input$groups_state
            sel<-input$groups_cell_clicked
            sel<-input$S1
            sel<<-sel
            # print(sel)
            sel
        # }
})

    # output$x1 = DT::renderDataTable(cars, server = FALSE)
    # 
    # # highlight selected rows in the scatterplot
    # output$x2 = renderPlot({
    #     s = input$x1_rows_selected
    #     par(mar = c(4, 4, 1, .1))
    #     plot(cars)
    #     if (length(s)) points(cars[s, , drop = FALSE], pch = 19, cex = 2)
    # })
    # 
    # # server-side processing
    # mtcars2 = mtcars[, 1:8]
    # output$x3 = DT::renderDataTable(mtcars2, server = TRUE)
    # 
    # print the selected indices
    output$examples_selected = renderPrint({
        s = input$examples_rows_selected
        if (length(s)) {
            cat('These rows were selected:\n\n')
            cat(s, sep = ', ')
        }
    })
    
    output$session_info<-renderPrint(sessionInfo())
    
}
