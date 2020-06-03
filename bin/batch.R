#!/usr/bin/R --no-save
wd<-sub("/bin$","",getwd())
setwd(wd)
source("shiny/global.R")
#FilesIn; GroupsSel; Exp; Rep    
tsize<-"2000"
Rep<-2
filesIn<-cbind(rf=FilesIn[,"file"],gr=unlist(GroupsSel[FilesIn[,"file"]]),wd="",name="",wf="",type="",ft="",ft1="",ft2="")
#test fastq?, convert to fastq
for(i in 1:nrow(filesIn)){
	rf <- filesIn[i,"rf"]
	ext<- tolower(sub('^.*[.$]',".",rf))
	s<-sub(ext,"",basename(rf))
	filesIn[i,"name"]<-s
	d<-paste0("data/",Exp,"/",s,"/")
	dir.create(paste0(d,"faTab"),recursive = T)
	filesIn[i,"wd"]<-d
	if(ext %in% c(".fq",".fastq")){
		filesIn[i,"wf"]<-rf
		filesIn[i,"type"]<-"fq"
		system(paste0("sed -n '1~4s/^@/>/p;2~4p' ",rf," | $HOME/conda/bin/fasta_formatter -t -o ",d,"faTab/",s,".faTab"),intern = TRUE)
	}
	if(ext %in% c(".fa",".fasta")){
		filesIn[i,"wf"]<-rf
		filesIn[i,"type"]<-"fa"
		system(paste0("$HOME/conda/bin/fasta_formatter -i ",rf," -t -o ",d,"faTab/",s,".faTab"),intern = TRUE)
	}
	if(ext %in% c(".bam",".sam",".cram")){
		system(paste0("samtools fastq ",rf," ",d,s,".fq" ))
		filesIn[i,"wf"]<-paste0(d,s,".fq")
		filesIn[i,"type"]<-"fq"
		system(paste0("sed -n '1~4s/^@/>/p;2~4p' ",d,s,".fq | $HOME/conda/bin/fasta_formatter -t -o ",d,"faTab/",s,".faTab"),intern = TRUE)
	}
	filesIn[i,"ft"]<-paste0(d,"faTab/",s,".faTab")
	for(r in 1:as.numeric(Rep)){
		system(paste0("shuf ",filesIn[i,"ft"]," | head -n ",tsize," > ",d,"faTab/",s,"_random",tsize,".",r,".faTab"))
		filesIn[i,paste0("ft",r)]<-paste0(d,"faTab/",s,"_random",tsize,".",r,".faTab")
	}
}

#blast per sample
#args <- c("R","--no-save",paste0("data/test"),paste0("../example-samples/Cancer_1.fa"),4)
for(i in 1:nrow(filesIn)){
	for(r in 1:as.numeric(Rep)){
		args <- c("sRNAflow","--no-save",filesIn[i,"wd"],paste0(filesIn[i,"name"],"_random",tsize,".",r),filesIn[i,paste0("ft",r)],4)
		source("bin/sRNAflow_blast_per_sample.R")
	}
}

#aggregate species  ########
source("bin/sRNAflow_aggregate_species.R")

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

