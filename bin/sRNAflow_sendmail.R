#!/usr/bin/R --no-save
library(sendmailR)
from <- sprintf("sRNAflow@biomed.lu.lv","sRNAflow")
to <- sprintf(email)
subject <- paste("sRNAflow",Exp)
body <- paste("sRNAflow",Exp)
bodyWithAttachment <- list(body)
zip(file.path(ED,paste0(Exp,"_fastQC.zip")),files=dir(file.path(ED,"qc"),".html",full.names = T),extras="-o -j -9")
zip(file.path(ED,paste0(Exp,"_isomiR-SEA.zip")),files=dir(ED,"_isomiR-SEA.xlsx",full.names = T,recursive = TRUE),extras="-o -j -9")
#TODO Limit size of mail to 20MB  ####
for(f in c(paste0(Exp,"_results.xlsx"),paste0(Exp,"_fastQC.zip"),paste0(Exp,"_isomiR-SEA.zip")))
	bodyWithAttachment<-append(bodyWithAttachment,mime_part(x=file.path(ED,f),name=f))
sendmail(from,to,subject,bodyWithAttachment,control=list(smtpServer=smtpServer))
