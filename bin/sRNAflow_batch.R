#!/usr/bin/env Rscript --vanilla 
# args <- commandArgs(TRUE)
wd<-sub("/bin$","",getwd())
setwd(wd)
source("shiny/global.R")

ED<-file.path(wd,"www","results",Exp)
dir.create(ED,recursive = TRUE, mode = "0777")
save(FilesIn,GroupsSel,Exp,specie,tsize,Rep,blast,qc,ad3,ad5,sizerange,lim,limS,log2FoldChange,padj,email,smtpServer,strategy,file = file.path(ED,"settings.RData"))
load(file.path(ED,"settings.RData"))
# FilesIn<-serverFiles[grep("-2.fq.gz",serverFiles[,"file"]),]; GroupsSel<-rep("ignore",nrow(FilesIn)); names(GroupsSel)<-FilesIn[,"file"]

#download main genomes from Ensembl  ########
source("bin/sRNAflow_downloadMainGenomes.R")

#Data trimming and QC  ########
source("bin/sRNAflow_filesIn_subsets.R")
load(file.path(ED,"filesIn.RData"))

#BLAST analysis  ########
source(file.path(wd,"bin","sRNAflow_blast_per_sample.R"))

#aggregate species  ########
source("bin/sRNAflow_aggregate_species.R")

#download genomes and index generation   #####
source("bin/sRNAflow_downloadGenomes.R")

#Mapping & RNA types catalog  ####
err<-foreach(i=1:nrow(filesIn)) %dopar% { system(paste(file.path(wd,"bin","sRNAflow_mapping_and_RNA_catalog.sh"),
	"-s",strategy,"-v",specie,"-n",filesIn[i,"name"],"-r",filesIn[i,"rf"],"-f",filesIn[i,"wf"],"-t",filesIn[i,"type"],"-o",ED),intern = TRUE); }

#Differential expression analysis & report preparation ########
source("bin/sRNAflow_DESeq2.R")

#Genomeless analysis  ########
source("bin/sRNAflow_genomeless.R")

#Send results  ########
source("bin/sRNAflow_sendmail.R")
