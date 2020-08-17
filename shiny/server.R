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
        if(!dir.exists(file.path(wd,"www","upload"))) dir.create(file.path(wd,"www","upload"),recursive = TRUE)
        out<-gsub(" ","_",input$file_upload$name)
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
        observeEvent(input$refresh_examples,{
            examples<-dir(path=examplesdir,pattern = ".(fastq|fq|fasta|fa|bam|sam|cram)(|.(gz|bz2|xz))$", full.names = FALSE, recursive = TRUE, include.dirs = TRUE)
            if(length(examples)>0){
                examples<-cbind(file=file.path("example-samples",examples),
                                size=humanReadable(file.info(file.path(examplesdir,examples))$size),
                                date=format(file.info(file.path(examplesdir,examples))$mtime,"%d.%m.%Y %H:%M:%OS"))
            } else examples<-rbind(rep(NA,3))[-1,]
        })
        examples<<-examples
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
        save(FilesIn,file=file.path(wd,"www","db","FilesIn.RData"))
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

    output$sel = renderTable({
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
            save(GroupsSel,file=file.path(wd,"www","db","GroupsSel.RData"))
        }
        # table(unlist(sel))
        rbind(number=table(Number=as.character(GroupsSel) )," "=c(" ") )
    })
    output$session_info<-renderPrint(sessionInfo())
    output$Config<-renderText({
        Exp       <<-gsub(" ","_",trimws(input$Exp))
        specie    <<-input$specie
        tsize     <<-input$tsize
        Rep       <<-input$Rep
        blast     <<-input$blast
        qc        <<-input$qc
        ad3       <<-input$ad3
        ad5       <<-input$ad5
        sizerange <<-input$sizerange
        lim       <<-input$lim
        log2FoldChange<<-input$log2FoldChange
        padj      <<-input$padj
        email     <<-input$email
        smtpServer<<-input$smtpServer
        save(Exp,specie,tsize,Rep,blast,qc,ad3,ad5,sizerange,lim,log2FoldChange,padj,email,smtpServer,file=file.path(wd,"www","db","Config.RData"))
        # list(wd=wd,getwd=getwd(),Exp=Exp,specie=specie,tsize=tsize,Rep=Rep,blast=blast,ad3=ad3,ad5=ad5,sizerange =sizerange,lim =lim,log2FoldChange=log2FoldChange,padj =padj,email =email,smtpServer=smtpServer)
        Exp
    })
    
    observeEvent(input$start,{
        ED<-file.path(wd,"www","results",Exp)
        if(dir.exists(ED)){
            showModal(modalDialog(
                title = "Warning",
                "Analysis with this name already started!"
            ))
        } else {
            # parallel::mcfork(source("bin/batch.R"),estranged = TRUE)
            # parallel::mcparallel(source("bin/batch.R"),detached = TRUE)
            # parallel::mccollect(wait = F)
            p <- parallel:::mcfork(estranged = TRUE)
            setwd(wd)
            if (inherits(p, "masterProcess")) {
                source("bin/batch.R")
                parallel:::mcexit()
            }
            withProgress(message = 'Starting job', {
                for(i in 1:4){
                    Sys.sleep(1)
                    incProgress(1/4)
                }
            })
        }
    })
    observeEvent(input$download_examples,{
        p <- parallel:::mcfork(estranged = TRUE)
        setwd(wd)
        if (inherits(p, "masterProcess")) {
            source("bin/sRNAflow_example_set.R")
            parallel:::mcexit()
        }
    })
    #         # system(paste("Rscript","bin/batch.R",Exp),wait = FALSE)
    #         system(paste0("R -e 'Exp<-\"",Exp,"\"; source(\"bin/batch.R\")'"),wait = FALSE)
    #         # source("bin/batch.R")
}




