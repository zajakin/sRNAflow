#!/usr/bin/R --no-save
wd<-sub("/bin$","",getwd())
setwd(wd)
source("shiny/global.R")
FilesIn
GroupsSel
args <- c("R","--no-save",paste0("data/test"),paste0("../example-samples/Cancer_1.fa"),4)
source("bin/parse_blast_output.R")


for(file in FilesIn[,"file"]) system(paste("bin/mapping_and_RNA_catalog.sh",file),intern = TRUE)


DESeq2

GOstats

miRanda

MultiQC

mail

