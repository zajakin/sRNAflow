#!/usr/bin/R --no-save
wd<-sub("/bin$","",getwd())
setwd(wd)
source("shiny/global.R")
#FilesIn; GroupsSel; Exp;
filesIn<-cbind(rf=FilesIn[,"file"],gr=unlist(GroupsSel[FilesIn[,"file"]]),fq="",fa="")
#test fastq?, convert to fastq
for(i in 1:nrow(filesIn)){
	rf <- filesIn[i,"rf"]
	ext<- tolower(sub('^.*[.$]',".",rf))
	s<-sub(ext,"",basename(rf))
	d<-paste0("data/",Exp,"/",s,"/")
	dir.create(d,recursive = T)
	if(ext %in% c(".fq",".fastq")) filesIn[i,"fq"]<-rf
	if(ext %in% c(".fa",".fasta")) filesIn[i,"fa"]<-rf
	if(ext %in% c(".bam",".sam",".cram")){
		system(paste0("samtools fastq ",rf," ",d,s,".fq" ))
		filesIn[i,"fq"]<-paste0(d,s,".fq" )
	}
}

#blast per sample
#args <- c("R","--no-save",paste0("data/test"),paste0("../example-samples/Cancer_1.fa"),4)
args <- c("sRNAflow","--no-save",paste0("data/",Exp),FilesIn[1,1],4)
source("bin/parse_blast_output.R")

#aggregate species

#download genomes

#index generation

#mapping
for(file in FilesIn[,"file"]) system(paste("bin/mapping_and_RNA_catalog.sh",file),intern = TRUE)

#RNA catalog

#Diff expression DESeq2

#GOstats

#miRanda

#MultiQC

#mail

