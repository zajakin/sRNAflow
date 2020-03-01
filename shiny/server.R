# autoInvalidate <- reactiveTimer(10000)
# observe({
#     autoInvalidate()
#     cat(".")
# })

server <- function(input, output) {
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
    # output$clock <- renderText({
    #     invalidateLater(1000)
    #     format(Sys.time(),"%d.%m.%Y %H:%M:%OS")
    # })
    output$filesUpload <- renderText({
        req(input$file_upload)
        if(!dir.exists(paste0(wd,"/data/input"))) dir.create(paste0(wd,"/data/input"),recursive = TRUE)
        out<-gsub(" ","_",paste0(wd,"/data/input/",input$file_upload$name))
        file.copy(input$file_upload$datapath,out)
        file.remove(input$file_upload$datapath)
        out
    })

    output$filesUploaded = DT::renderDataTable({
        if (!is.null(input$file_upload)) filesUploaded<<-rbind(filesUploaded,filesUploadRow(input$file_upload))
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
    
    output$filesIn = DT::renderDataTable({
        filesIn <- FilesIn
        s = input$filesIn_rows_selected
        if (length(s)) filesIn <- rbind(filesIn[-s,])
        s = input$filesUploaded_rows_selected
        observeEvent(input$filesUploaded_rows_selected,shinyjs::reset("filesUploaded_rows_selected"))
        if (length(s)) filesIn <- rbind(filesUploaded[s,],filesIn)
        s = input$examples_rows_selected
        observeEvent(input$examples_rows_selected,shinyjs::reset("examples_rows_selected"))
        if (length(s)) filesIn <- rbind(examples[s,],filesIn)
        s = input$serverFiles_rows_selected
        observeEvent(input$serverFiles_rows_selected,shinyjs::reset("serverFiles_rows_selected"))
        if (length(s)) filesIn <- rbind(serverFiles[s,],filesIn)
        if(nrow(filesIn)>0) colnames(filesIn)<-c("file","size","date")
        filesIn <- unique(filesIn)
        FilesIn <<- filesIn
        save(FilesIn,file=paste0(wd,"/data/FilesIn.RData"))
        filesIn
    }, server = TRUE)

    output$groups = DT::renderDataTable({
        tmp<-length(c(input$filesUploaded_rows_selected,input$examples_rows_selected,input$serverFiles_rows_selected,input$filesIn_rows_selected))>0
        groups <- cbind(rbind(FilesIn),test="",control="",ignore="")
        colnames(groups)<-c("file","size","date","test","control","ignore")
        if(nrow(groups)>0){
            # rownames(groups)<-paste0("S",1:nrow(groups))
            groups[,"test"]<-sprintf('<input type="radio" name="%s" value="1"/>',groups[,"file"])
            groups[,"control"]<-paste0('<input type="radio" name="',groups[,"file"],'" value="0"/>')
            groups[,"ignore"] <-paste0('<input type="radio" name="',groups[,"file"],'" value="-1"/>')
        }
        Groups <<- groups
        groups
    }, server = FALSE, escape = FALSE, selection = 'none', options = list(dom = 't', paging = FALSE, ordering = FALSE), 
    callback = DT::JS("table.rows().every(function(i, tab, row) {
          var $this = $(this.node());
          $this.attr('id', this.data()[0]);
          $this.addClass('shiny-input-radiogroup');
        });
        Shiny.unbindAll(table.table().node());
        Shiny.bindAll(table.table().node());")
    )

    output$sel = renderPrint({
        if(length(input$groups_cell_clicked)>0){
            # print(input$groups_cell_clicked)
            sel<-c()
            sel<-sapply(Groups[,"file"], function(i) c<-c(sel,input[[i]]))
            names(sel)<-Groups[,"file"]
            Sel <<- sel
            save(Sel,file=paste0(wd,"/data/GroupsSel.RData"))
            sel
        }
    })

    output$examples_selected = renderPrint({
        s = input$examples_rows_selected
        if (length(s)) {
            cat('These rows were selected:\n\n')
            cat(s, sep = ', ')
        }
    })
    
    output$session_info<-renderPrint(sessionInfo())
}
