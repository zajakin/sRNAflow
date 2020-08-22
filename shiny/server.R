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
        if(!dir.exists(file.path(wd,"www","upload"))) dir.create(file.path(wd,"www","upload"),recursive = TRUE, mode = "0777")
        out<-gsub(" ","_",input$file_upload$name)
        file.copy(input$file_upload$datapath,file.path(wd,"www","upload",out))
        file.remove(input$file_upload$datapath)
        out
    })

    output$filesUploaded = DT::renderDataTable({
        if (!is.null(input$file_upload)) filesUploaded<<-rbind(filesUploaded,filesUploadRow(input$file_upload))
        filesUploaded
    }, server = TRUE)

    output$examples = DT::renderDataTable({
        examplesdir<-file.path(wd,"www","upload","example-samples")
        if(!dir.exists(examplesdir)) dir.create(examplesdir,recursive = TRUE, mode = "0777")
        examples<-dir(path=examplesdir,pattern = ".(fastq|fq|fasta|fa|bam|sam|cram)(|.(gz|bz2|xz))$", full.names = FALSE, recursive = TRUE, include.dirs = TRUE)
        if(length(examples)>0){
            examples<-cbind(file=file.path("example-samples",examples),
                            size=humanReadable(file.info(file.path(examplesdir,examples))$size),
                            date=format(file.info(file.path(examplesdir,examples))$mtime,"%d.%m.%Y %H:%M:%OS"))
        } else examples<-rbind(rep(NA,3))[-1,]
        if (length(input$filesIn_rows_selected)) observeEvent(input$examples_rows_selected,shinyjs::reset("examples_rows_selected"))
        input$refresh_examples
        colnames(examples) <- c("file","size","date")
        examples<<-examples
        examples
    }, server = TRUE)
    
    output$serverFiles = DT::renderDataTable({
        if(!dir.exists(file.path(wd,"www","upload"))) dir.create(file.path(wd,"www","upload"),recursive = TRUE, mode = "0777")
        serverFiles<-dir(path = file.path("www","upload"),pattern = ".(fastq|fq|fasta|fa|bam|sam|cram)(|.(gz|bz2|xz))$",full.names = FALSE, recursive = TRUE, include.dirs = TRUE)
        if(length(grep("^example-samples",serverFiles))>0) serverFiles<-serverFiles[-grep("^example-samples",serverFiles)]
        if(length(serverFiles)>0){
            serverFiles<-cbind(file=serverFiles,
                               size=humanReadable(file.info(file.path("www","upload",serverFiles))$size),
                               date=format(file.info(file.path("www","upload",serverFiles))$mtime,"%d.%m.%Y %H:%M:%OS"))
        } else serverFiles<-rbind(rep(NA,3))[-1,]
        if (length(input$filesIn_rows_selected)) observeEvent(input$serverFiles_rows_selected,shinyjs::reset("serverFiles_rows_selected"))
        length(input$refresh_serverFiles)
        colnames(serverFiles) <- c("file","size","date")
        serverFiles<<-serverFiles
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
        if(input$clear_filesIn){
            if(!exists("clear_filesIn_counter") || clear_filesIn_counter<input$clear_filesIn){ 
                filesIn <- rbind(rep(NA,3))[-1,]
                clear_filesIn_counter<<-input$clear_filesIn
                shinyjs::reset("clear_filesIn")
            }
        }
        if(nrow(filesIn)>0) colnames(filesIn)<-c("file","size","date")
        filesIn <- unique(filesIn)
        FilesIn <<- filesIn
        save(FilesIn,file=file.path(wd,"www","db","FilesIn.RData"))
        filesIn
    }, server = TRUE)

    output$groups = DT::renderDataTable({
        tmp<-length(c(input$filesUploaded_rows_selected,input$examples_rows_selected,input$serverFiles_rows_selected,input$filesIn_rows_selected))>0
        if(nrow(rbind(FilesIn))>0){
            groups <- cbind(rbind(FilesIn),test=c(""),control=c(""),ignore=c(""))
        } else groups <- rbind(rep(0,6))[-1,]
        choices <- c("test","control","ignore")
        colnames(groups)<-c("file","size","date",choices)
        if(nrow(groups)>0) groups[,4:6]<-paste0('<input type="radio" name="',groups[,"file"],'" value="',rep(choices,each=nrow(groups)),'"/>')
        if(nrow(rbind(FilesIn))>0){
            GroupsSel<<-GroupsSel[names(GroupsSel) %in% FilesIn[,"file"]]
            for(s in names(GroupsSel)) if(!is.null(GroupsSel[s][[1]])) groups[groups[,"file"]==s,GroupsSel[s][[1]]]<-sub("/>"," checked/>",groups[groups[,"file"]==s,GroupsSel[s][[1]]])
        }
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

    output$reports = DT::renderDataTable({
        input$refresh_reports
        exps<-dir(file.path(wd,"www","results"))
        reports <- cbind(xlsx=c(""),fastQC=c(""),multiQC=c(""),"isomiR-SEA"=c(""))[-1,]
        for(i in exps){
            row<-c()
            for(j in c("_results.xlsx","_fastQC.zip","_multiqc.html","_isomiR-SEA.zip")) 
                if(file.exists(file.path(wd,"www","results",i,paste0(i,j)))){
                    # row<-c(row,a(href=paste0('/results/',i,'/',paste0(i,j)),paste0(i,j), download=NA, target="_blank"))
                    row<-c(row,paste0('<a href="/results/',i,'/',paste0(i,j),'" target="_blank">',paste0("...",j),'</a>'))
                } else row<-c(row,'<h2>-</h2>')
            reports<-rbind(reports,row)
            rownames(reports)[nrow(reports)]<-i
        }
        reports
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
            showModal(modalDialog(title = "Warning","Analysis with this name already started!"))
        } else {
            setwd(wd)
            system(paste0('Rscript --vanilla bin/batch.R "',wd,'"'), wait = FALSE)
            showModal(modalDialog('Analysis started. Go to "Reports" and press ""Refresh reports list"" to see generated files'))
            # if (inherits(p, "masterProcess")) {
            #     source("bin/batch.R")
            #     parallel:::mcexit()
            # }
            withProgress(message = 'Starting job', {
                for(i in 1:4){
                    Sys.sleep(1)
                    incProgress(1/4)
                }
            })
        }
    })
    observeEvent(input$download_examples,{
        setwd(wd)
        system(paste0('Rscript --vanilla bin/sRNAflow_example_set.R "',wd,'" "UF_"'), wait = FALSE)
        showModal(modalDialog('Downloading of 39 EV examples started.
                              Press "Refresh examples list" to see downloaded files'))
    })

    observeEvent(input$download_examples2,{
        setwd(wd)
        system(paste0('Rscript --vanilla bin/sRNAflow_example_set.R "',wd,'" "cells_"'), wait = FALSE)
        showModal(modalDialog('Downloading of 40 cells examples started.
                              Press "Refresh examples list" to see downloaded files'))
    })

    output$blastdb<-renderText({
        blastdb<-"No local BLAST database. Press 'Update local BLAST db' button to upload."
        if(file.exists(file.path(wd,"www","db","blast","db.done"))) blastdb<-"Local BLAST database exist and will be used in analysis."
        blastdb
    })

    observeEvent(input$blastdb,{
        setwd(wd)
        system(paste0('Rscript --vanilla bin/sRNAflow_local_blast_db.R "',wd,'"'), wait = FALSE)
        showModal(modalDialog("Downloading of local BLAST database started"))
    })

    observeEvent(input$gtfdb,{
        setwd(wd)
        system(paste0('Rscript --vanilla bin/sRNAflow_gtf_annotations.R "',wd,'"'), wait = FALSE)
        showModal(modalDialog("Downloading and merging of small RNA annotations started"))
    })
    #         # system(paste("Rscript","bin/batch.R",Exp),wait = FALSE)
    #         system(paste0("R -e 'Exp<-\"",Exp,"\"; source(\"bin/batch.R\")'"),wait = FALSE)
    #         # source("bin/batch.R")
}




