#!/usr/bin/R --no-save
filesIn<-cbind(rf=FilesIn[,"file"],gr=unlist(GroupsSel[FilesIn[,"file"]]),wd="",name="",wf="",type="",ft="",ft1="",ft2="")

trimm<-function(rf,ext,s,d,ad3,ad5,sizerange,arx){
	if(trimws(sub("#.*","",ad3)) != "") ad3<-paste("-a",trimws(sub("#.*","",ad3)))
	if(trimws(sub("#.*","",ad5)) != "") ad5<-paste("-g",trimws(sub("#.*","",ad5)))
	cat <- "cat "
	if(arx!="") cat <-"pigz -cd "
	qc<-""
	mark="^>"
	if(ext == "fastq"){ 
		qc<-"--quality-cutoff=20,20"
		mark="^@"
	}
	con<-file(paste0(d,"trimm.txt"),"wt")
	writeLines(paste0("File\t",s),con)
	writeLines(paste0("Raw\t",system(paste0(cat,rf," | grep -c \"",mark,"\" "),intern = TRUE)),con)

	system(paste("cutadapt --cores=0 --quiet",qc,ad3,ad5,"-m 1 -o ",paste0(d,s,"_r3.",ext),rf),intern = TRUE)
	writeLines(paste0("QC_and_adapter3\t",system(paste0("grep -c \"",mark,"\" ",paste0(d,s,"_r3.",ext)),intern = TRUE)),con)

	system(paste("cutadapt --cores=0 --quiet -M",sizerange[2],"-o",paste0(d,s,"_r4.",ext),paste0(d,s,"_r3.",ext)),intern = TRUE)
	writeLines(paste0("Filter_long\t",system(paste0("grep -c \"",mark,"\" ",paste0(d,s,"_r4.",ext)),intern = TRUE)),con)

	wf<-paste0(d,s,".",ext)
	system(paste("cutadapt --cores=0 --quiet -m",sizerange[1],"-o",wf,paste0(d,s,"_r4.",ext)),intern = TRUE)
	# writeLines(paste0("Filter_short\t",system(paste0("grep -c \"",mark,"\" ",wf),intern = TRUE)),con)
	if(ext == "fastq"){
		system(paste0("fastqc -q --threads ",core," -o ",ED,"/qc ",wf),intern = TRUE)
		system(paste0("sed -n '1~4s/^@/>/p;2~4p' ",wf," | $HOME/conda/bin/fasta_formatter -t -o ",d,'faTab/',s,'.faTab'),intern = TRUE)
	} else system(paste0("$HOME/conda/bin/fasta_formatter -i ",wf," -t -o ",d,'faTab/',s,'.faTab'))

	mstat<-'awk \'{ i++; seq[i] = length($NF); sum += length($NF) } END { asort(seq,sorted,"@val_num_asc"); print("Filter_short\\t" i "\\nMedian_reads_length\\t" sorted[int(i/2)] "\\nMean_reads_length\\t" sum / i) }\' '
	writeLines(system(paste0(mstat,d,'faTab/',s,'.faTab'),intern=TRUE),con)
	close(con)
	file.remove(paste0(d,s,"_r3.",ext),paste0(d,s,"_r4.",ext))
	system(paste("pigz -9f ",wf))
	wf<-paste0(wf,".gz")
	return(wf)
}

filesIn<-foreach(i=1:nrow(filesIn),.combine = rbind) %dopar% {
	rf <- filesIn[i,"rf"]
	ext<- tolower(sub('^.*[.$]',".",rf))
	arx<- ""
	if(ext %in% c(".gz",".bz2",".xz")){
		arx<- sub('^.*[.$]',".",rf)
		ext<- tolower(sub('^.*[.$]',".",sub(arx,'',rf)))
	}
	s<-sub(paste0(ext,arx),"",basename(rf))
	filesIn[i,"name"]<-s
	d<-file.path("www",Exp,s,"")
	if(!dir.exists(paste0(d,"faTab"))) dir.create(paste0(d,"faTab"),recursive = T)
	filesIn[i,"wd"]<-d
	filesIn[i,"ft"]<-paste0(d,"faTab/",s,".faTab")
	if(ext %in% c(".fa",".fasta")){
		filesIn[i,"wf"]<-trimm(rf,"fasta",s,d,ad3,ad5,sizerange,arx)
		filesIn[i,"type"]<-"fa"
		# system(paste0("pigz -cd ",filesIn[i,"wf"]," | $HOME/conda/bin/fasta_formatter -t -o ",filesIn[i,"ft"]),intern = TRUE)
	}
	if(ext %in% c(".bam",".sam",".cram")){
		system(paste0("samtools fastq --threads ",core," ",rf," -o ",d,s,".fq.gz" ))
		filesIn[i,"rf"]<-rf<-paste0(d,s,".fq.gz")
		arx<-".gz"
		ext<-".fq"
	}
	if(ext %in% c(".fq",".fastq")){
		filesIn[i,"wf"]<-trimm(rf,"fastq",s,d,ad3,ad5,sizerange,arx)
		filesIn[i,"type"]<-"fq"
		# system(paste0("pigz -cd ",filesIn[i,"wf"]," | sed -n '1~4s/^@/>/p;2~4p' | $HOME/conda/bin/fasta_formatter -t -o ",filesIn[i,"ft"]),intern = TRUE)
	}
	for(r in 1:as.numeric(Rep)){
		# system(paste0("shuf -n ",tsize," ",filesIn[i,"ft"]," -o ",d,"faTab/",s,"_random",tsize,".",r,".faTab"), intern = T)
		filesIn[i,paste0("ft",r)]<-paste0(d,"faTab/",s,"_random",tsize,".",r,".faTab")
	}
	file.remove(filesIn[i,"ft"])
	filesIn[i,]
}
save(filesIn,file = file.path(wd,"www",Exp,"filesIn.RData"))
