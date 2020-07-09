#!/usr/bin/R --no-save
filesIn<-cbind(rf=FilesIn[,"file"],gr=unlist(GroupsSel[FilesIn[,"file"]]),wd="",name="",wf="",type="",ft="",ft1="",ft2="")
#test fastq?, convert to fastq
for(i in 1:nrow(filesIn)){
	rf <- filesIn[i,"rf"]
	ext<- tolower(sub('^.*[.$]',".",rf))
	s<-sub(ext,"",basename(rf))
	filesIn[i,"name"]<-s
	d<-paste0("data/",Exp,"/",s,"/")
	if(!dir.exists(paste0(d,"faTab"))) dir.create(paste0(d,"faTab"),recursive = T)
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
#TODO Shuf show error
		filesIn[i,paste0("ft",r)]<-paste0(d,"faTab/",s,"_random",tsize,".",r,".faTab")
	}
}
save(filesIn,file = file.path(wd,"data",Exp,"filesIn.RData"))
