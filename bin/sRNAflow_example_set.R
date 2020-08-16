#!/usr/bin/R --no-save
WD<-file.path(wd,"www","upload","example-samples")
desc<-file.path(WD,"filereport_read_run_PRJNA293274.tsv")
dir.create(WD,recursive=TRUE)
download.file("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA293274&result=read_run&fields=fastq_ftp,sample_title&format=tsv&download=true",desc)
tab<-read.table(desc, sep = "\t",header = F,skip = 1,comment.char = "",fill = T)
colnames(tab)<-read.table(desc,header = FALSE,nrows = 1,skip = 0,comment.char = "")
tab<-tab[grep("UF_",tab[,"sample_title"]),]
for(i in 1:nrow(tab)) download.file(paste0("ftp://",tab[i,"fastq_ftp"]),file.path(WD,paste0(tab[i,"sample_title"],".fastq.gz")))
