library(shiny)
library(gdata)
options(shiny.maxRequestSize=10000*1024^2)

wd<-sub("/shiny$","",getwd())
# source(paste0(wd,"/shiny/utils.R"))
if(!dir.exists(paste0(wd,"/data/example-samples"))) dir.create(paste0(wd,"/data/example-samples"),recursive = TRUE)
examples<-dir(path = paste0(wd,"/data/example-samples"),full.names = FALSE, recursive = TRUE, include.dirs = TRUE)
if(length(examples)>0)
  examples<-cbind(file=paste0("data/example-samples/",examples),
				size=humanReadable(file.info(paste0(wd,"/data/example-samples/",examples))$size),
				date=format(file.info(paste0(wd,"/data/example-samples/",examples))$mtime,"%d.%m.%Y %H:%M:%OS"))
if(!dir.exists(paste0(wd,"/data/input"))) dir.create(paste0(wd,"/data/input"),recursive = TRUE)
serverFiles<-dir(path = paste0(wd,"/data/input"),full.names = FALSE, recursive = TRUE, include.dirs = TRUE)
if(length(serverFiles)>0)
  serverFiles<-cbind(file=paste0("data/input/",serverFiles),
				size=humanReadable(file.info(paste0(wd,"/data/input/",serverFiles))$size),
				date=format(file.info(paste0(wd,"/data/input/",serverFiles))$mtime,"%d.%m.%Y %H:%M:%OS"))
filesIn           <- filesUploaded           <- rbind(rep(NA,3))[-1,]
colnames(filesIn) <- colnames(filesUploaded) <- c("file","size","date")

# timeout<-2147483
# options(app_init_timeout=timeout)
# options(shiny.app_init_timeout=timeout)
# options(app_idle_timeout=timeout)
# options(shiny.app_idle_timeout=timeout)
# options(app_connection_timeout=timeout)
# options(shiny.app_connection_timeout=timeout)
# options(app_read_timeout=timeout)
# options(shiny.app_read_timeout=timeout)
# options(shiny.trace=T)
# options(port=4110)

