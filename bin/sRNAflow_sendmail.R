#!/usr/bin/R --no-save
sendresults<-function(){
library(sendmailR)
from <- sprintf("sRNAflow@biomed.lu.lv")
to <- sprintf(email)
subject <- paste("sRNAflow",Exp)
body <- paste("sRNAflow",Exp)
bodyWithAttachment <- list(body)
system(paste("cd ",ED,"; multiqc ."))
file.rename(file.path(ED,"multiqc_report.html"),file.path(ED,paste0(Exp,"_multiqc.html")))
mailsize<-1000
for(f in c(paste0(Exp,"_results.xlsx"),paste0(Exp,"_fastQC.zip"),paste0(Exp,"_multiqc.html"),paste0(Exp,"_isomiR-SEA.zip"))){
	if( (mailsize+file.size(file.path(ED,f)))>1.8e+7 ) next;
	bodyWithAttachment<-append(bodyWithAttachment,mime_part(x=file.path(ED,f),name=f))
}
sendmail(from,to,subject,bodyWithAttachment,control=list(smtpServer=smtpServer))
}

print(paste(date(),"Send results"))
sendresults()
