#!/usr/bin/R --no-save
# args <- commandArgs(TRUE)
wd<-sub("/bin$","",getwd())
setwd(wd)
source("shiny/global.R")

ED<-file.path(wd,"www","results",Exp)
dir.create(ED,recursive = TRUE, mode = "0777")
save(FilesIn,GroupsSel,Exp,specie,tsize,Rep,blast,qc,ad3,ad5,sizerange,lim,log2FoldChange,padj,email,smtpServer,strategy,file = file.path(ED,"settings.RData"))
load(file.path(ED,"settings.RData"))
# FilesIn<-serverFiles[grep("-2.fq.gz",serverFiles[,"file"]),]; GroupsSel<-rep("ignore",nrow(FilesIn)); names(GroupsSel)<-FilesIn[,"file"]

#download main genomes from Ensembl  ########
source("bin/sRNAflow_downloadMainGenomes.R")

source("bin/sRNAflow_filesIn_subsets.R")
load(file.path(ED,"filesIn.RData"))

source(file.path(wd,"bin","sRNAflow_blast_per_sample.R"))
comb<-cbind(rep(1:nrow(filesIn),each=as.numeric(Rep)),rep(1:as.numeric(Rep),times=nrow(filesIn)))
err<-foreach(combr=1:nrow(comb),.verbose = T) %dopar% blast_per_sample(idr=comb[combr,1],re=comb[combr,2],wd,filesIn,tsize,core,ED)

#aggregate species  ########
source("bin/sRNAflow_aggregate_species.R")

#download genomes and index generation   #####
source("bin/sRNAflow_downloadGenomes.R")

#Mapping & RNA types catalog  ####
tax<-"_tax_filtExactRC"
err<-foreach(i=1:nrow(filesIn)) %dopar% { system(paste(file.path(wd,"bin","sRNAflow_mapping_and_RNA_catalog.sh"),
	"-s",strategy,"-v",specie,"-n",filesIn[i,"name"],"-r",filesIn[i,"rf"],"-f",filesIn[i,"wf"],"-t",filesIn[i,"type"],"-o",ED),intern = TRUE); }

source("bin/sRNAflow_DESeq2.R")

source("bin/sRNAflow_sendmail.R")
