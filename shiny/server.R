# autoInvalidate <- reactiveTimer(10000)
# observe({
#     autoInvalidate()
#     cat(".")
# })

server <- function(input, output, session) {
    output$messageMenu <- renderMenu({
        # Code to generate each of the messageItems here, in a list. This assumes that messageData is a data frame with two columns, 'from' and 'message'.
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
        choices <- c("test","control","ignore")
        colnames(groups)<-c("file","size","date",choices)
        if(nrow(groups)>0) groups[,4:6]<-paste0('<input type="radio" name="',groups[,"file"],'" value="',rep(choices,each=nrow(groups)),'"/>')
        GroupsSel<<-GroupsSel[names(GroupsSel) %in% FilesIn[,"file"]]
        for(s in names(GroupsSel)) if(!is.null(GroupsSel[s][[1]])) groups[groups[,"file"]==s,GroupsSel[s][[1]]]<-sub("/>"," checked/>",groups[groups[,"file"]==s,GroupsSel[s][[1]]])
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
        sel <- GroupsSel
        if(length(c(input$groups_cell_clicked,input$filesUploaded_rows_selected,input$examples_rows_selected,input$serverFiles_rows_selected,input$filesIn_rows_selected))>0){
            # browser()
            # View(input$groups_cell_clicked)
            sel<-c()
            for(s in FilesIn[,"file"]){
                if(length(grep(s,names(input)))>0){
                    sel[s]<-as.character(input[[s]])
                    names(sel)[length(sel)]<-s
                }
            }
            sel<-c(GroupsSel[!names(GroupsSel) %in% names(sel)],sel)
            GroupsSel <<- sel
            save(GroupsSel,file=paste0(wd,"/data/GroupsSel.RData"))
        }
        sel
    })
    output$session_info<-renderPrint(sessionInfo())
    output$Config<-renderPrint({
        Exp<<-input$Exp
        save(Exp,file=paste0(wd,"/data/Config.RData"))
        Exp
    })
}
