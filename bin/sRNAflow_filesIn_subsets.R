#!/usr/bin/R --no-save
filesIn<-cbind(rf=FilesIn[,"file"],gr=unlist(GroupsSel[FilesIn[,"file"]]),wd="",name="",wf="",type="",ft="",ft1="",ft2="")

trimm<-function(rf,ext,s,d,ad3,ad5,sizerange){
	if(trimws(sub("#.*","",ad3)) != "") ad3<-paste("-a",trimws(sub("#.*","",ad3)))
	if(trimws(sub("#.*","",ad5)) != "") ad5<-paste("-g",trimws(sub("#.*","",ad5)))
	qc<-""
	if(ext == "fastq") qc<-"--quality-cutoff=20,20"
	
	# fasta_formatter -i $ff -t > $out/$f/${f}.faTab
	# gawk 'length($NF) < 43' $out/$f/${f}.faTab > $out/$f/${f}_r4.faTab
	# gawk 'length($NF) > 9 {print ">" $1 "\n" $NF}' $out/$f/${f}_r4.faTab > $out/$f/${f}.fasta

	# fastqc -o $out/qc_raw $ff  > /dev/null 2>&1
	system(paste("cutadapt --cores=0 --quiet",qc,ad3,ad5,"-m 1 -o ",paste0(d,s,"_r3.",ext),rf))
	system(paste("cutadapt --cores=0 --quiet -M",sizerange[2],"-o",paste0(d,s,"_r4.",ext),paste0(d,s,"_r3.",ext)))
	wf<-paste0(d,s,".",ext)
	system(paste("cutadapt --cores=0 --quiet -m",sizerange[1],"-o",wf,paste0(d,s,"_r4.",ext)))
	return(wf)
}

#test fastq?, convert to fastq
filesIn<-foreach(i=1:nrow(filesIn),.combine = rbind) %dopar% {
	rf <- filesIn[i,"rf"]
	ext<- tolower(sub('^.*[.$]',".",rf))
	s<-sub(ext,"",basename(rf))
	filesIn[i,"name"]<-s
	d<-file.path("www",Exp,s,"")
	if(!dir.exists(paste0(d,"faTab"))) dir.create(paste0(d,"faTab"),recursive = T)
	filesIn[i,"wd"]<-d
	if(ext %in% c(".fa",".fasta")){
		filesIn[i,"wf"]<-trimm(rf,"fasta",s,d,ad3,ad5,sizerange)
		filesIn[i,"type"]<-"fa"
		system(paste0("$HOME/conda/bin/fasta_formatter -i ",filesIn[i,"wf"]," -t -o ",d,"faTab/",s,".faTab"),intern = TRUE)
	}
	if(ext %in% c(".bam",".sam",".cram")){
		system(paste0("samtools fastq ",rf," ",d,s,".fq" ))
		filesIn[i,"rf"]<-rf<-paste0(d,s,".fq")
		ext<-".fq"
	}
	if(ext %in% c(".fq",".fastq")){
		filesIn[i,"wf"]<-trimm(rf,"fastq",s,d,ad3,ad5,sizerange)
		filesIn[i,"type"]<-"fq"
		system(paste0("sed -n '1~4s/^@/>/p;2~4p' ",filesIn[i,"wf"]," | $HOME/conda/bin/fasta_formatter -t -o ",d,"faTab/",s,".faTab"),intern = TRUE)
	}
	filesIn[i,"ft"]<-paste0(d,"faTab/",s,".faTab")
	for(r in 1:as.numeric(Rep)){
		system(paste0("shuf -n ",tsize," ",filesIn[i,"ft"]," -o ",d,"faTab/",s,"_random",tsize,".",r,".faTab"), intern = T)
		filesIn[i,paste0("ft",r)]<-paste0(d,"faTab/",s,"_random",tsize,".",r,".faTab")
	}
	filesIn[i,]
}
save(filesIn,file = file.path(wd,"www",Exp,"filesIn.RData"))
