#!/usr/bin/R --no-save
filesIn<-cbind(rf=FilesIn[,"file"],gr=unlist(GroupsSel[FilesIn[,"file"]]),wd="",name="",wf="",type="",ft="")
for(r in 1:as.numeric(Rep)){ filesIn<-cbind(filesIn,ftN=""); colnames(filesIn)[ncol(filesIn)]<-paste0("ft",r); }
if(!dir.exists(file.path(ED,"qc"))) dir.create(file.path(ED,"qc"),recursive = T, mode = "0777")
if(!dir.exists(file.path(ED,"species_diagrams"))) dir.create(file.path(ED,"species_diagrams"),recursive = TRUE, mode = "0777")

trimm<-function(rf,ext,s,d,qc,ad3,ad5,sizerange,arx,wf=rf){
	ad3<-trimws(sub("#.*","",ad3))
	ad5<-trimws(sub("#.*","",ad5))
	if(ad3 != "") ad3<-paste("-a",ad3)
	if(ad5 != "") ad5<-paste("-g",ad5)
	cat <- "cat "
	if(arx!="") cat <-"pigz -cd "
	qcc<-""
	mark="^>"
	if(ext == "fastq"){ 
		qcc<-paste0("--quality-cutoff=",qc,",",qc)
		mark="^@"
	}
	con<-file(file.path(d,"logs","trimm.txt"),"wt")
	writeLines(paste0("File\t",s),con)
	writeLines(paste0("Raw\t",system(paste0(cat,rf," | grep -c \"",mark,"\" "),intern = TRUE)),con)

	system(paste("cutadapt --cores=0 --quiet",qcc,ad3,ad5,"-m 1 -o ",paste0(d,s,"_r3.",ext),rf),intern = TRUE)
	writeLines(paste0("QC_and_adapter3\t",system(paste0("grep -c \"",mark,"\" ",paste0(d,s,"_r3.",ext)),intern = TRUE)),con)

	system(paste("cutadapt --cores=0 --quiet -M",sizerange[2],"-o",paste0(d,s,"_r4.",ext),paste0(d,s,"_r3.",ext)),intern = TRUE)
	writeLines(paste0("Filter_long\t",system(paste0("grep -c \"",mark,"\" ",paste0(d,s,"_r4.",ext)),intern = TRUE)),con)


	system(paste0("bash -c \"sed -n '1~4N;2~4N;3~4N;s/\\n/	/g;4~4p' ",paste0(d,s,"_r4.",ext)," | grep -v -F -f ",ED,"/exclude | tr '\t' '\n' \" > ",paste0(d,s,"_r5.",ext)),intern = TRUE)
	writeLines(paste0("Environment\t",system(paste0("grep -c \"",mark,"\" ",paste0(d,s,"_r5.",ext)),intern = TRUE)),con)

	wf<-paste0(d,s,".",ext)
	system(paste("cutadapt --cores=0 --quiet -m",sizerange[1],"-o",wf,paste0(d,s,"_r5.",ext)),intern = TRUE)

	fasta2tab<-" | gawk '/^>/ {printf(\"%s%s\\t\",(N>0?\"\\n\":\"\"),$1);N++;next;} {printf(\"%s\",$0);} END {printf(\"\\n\");}' > "
	if(ext == "fastq"){
		system(paste0("fastqc -q --threads ",core," -o ",ED,"/qc ",wf),intern = TRUE)
		system(paste0("sed -n '1~4s/^@/>/p;2~4p' ",wf,fasta2tab,d,'faTab/',s,'.faTab'),intern = TRUE)
	} else system(paste0("cat ",wf,fasta2tab,d,'faTab/',s,'.faTab'))
	system(paste0("bowtie2 --time --end-to-end -k 201 -p ",core," --mm -x ",file.path(wd,"www","db","genomes","bowtie2","univec","univec"),
		  " --un /dev/null --no-unal -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -U ",wf," -S ",d,"univec.sam > ",d,"logs/univec.log 2>&1"))
	writeLines(paste0("UniVec\t",system(paste0("grep \"overall alignment rate\" ",d,"logs/univec.log | gawk '{print $1}'"),intern = TRUE)),con)
	
	mstat<-'gawk \'{ i++; seq[i] = length($NF); sum += length($NF) } END { asort(seq,sorted,"@val_num_asc"); print("Filter_short\\t" i "\\nMedian_reads_length\\t" sorted[int(i/2)] "\\nMean_reads_length\\t" sum / i) }\' '
	writeLines(system(paste0(mstat,d,'faTab/',s,'.faTab'),intern=TRUE),con)
	close(con)
	file.remove(paste0(d,s,"_r3.",ext),paste0(d,s,"_r4.",ext),paste0(d,s,"_r5.",ext))
	system(paste("pigz -9f ",wf))
	wf<-paste0(wf,".gz")
	return(wf)
}

filesInGen<-function(i){
	rf <- file.path(wd,"www","upload",filesIn[i,"rf"])
	ext<- tolower(sub('^.*[.$]',".",rf))
	arx<- ""
	if(ext %in% c(".gz",".bz2",".xz")){
		arx<- sub('^.*[.$]',".",rf)
		ext<- tolower(sub('^.*[.$]',".",sub(arx,'',rf)))
	}
	s<-sub(paste0(ext,arx),"",basename(rf))
	filesIn[i,"name"]<-s
	d<-file.path(ED,s,"")
	if(!dir.exists(paste0(d,"logs"))) dir.create(paste0(d,"logs"),recursive = T)
	if(!dir.exists(paste0(d,"faTab"))) dir.create(paste0(d,"faTab"),recursive = T)
	filesIn[i,"wd"]<-d
	filesIn[i,"ft"]<-paste0(d,"faTab/",s,".faTab")
	if(ext %in% c(".fa",".fasta")){
		filesIn[i,"wf"]<-trimm(rf,"fasta",s,d,qc,ad3,ad5,sizerange,arx)
		filesIn[i,"type"]<-"fa"
	}
	if(ext %in% c(".bam",".sam",".cram")){
		system(paste0("samtools fastq --threads ",core," ",rf," -o ",d,s,".fq.gz" ))
		filesIn[i,"rf"]<-rf<-paste0(d,s,".fq.gz")
		arx<-".gz"
		ext<-".fq"
	}
	if(ext %in% c(".fq",".fastq")){
		filesIn[i,"wf"]<-trimm(rf,"fastq",s,d,qc,ad3,ad5,sizerange,arx)
		filesIn[i,"type"]<-"fq"
	}
	for(r in 1:as.numeric(Rep)){
		system(paste0("shuf -n ",tsize," ",filesIn[i,"ft"]," -o ",d,"faTab/",s,"_random",tsize,".",r,".faTab"), intern = T)
		filesIn[i,paste0("ft",r)]<-paste0(d,"faTab/",s,"_random",tsize,".",r,".faTab")
	}
	system(paste0("cut -f 2 ",filesIn[i,"ft"]," | sort | uniq -c | sort -r | pigz -9cf  > ",d,s,".reads.gz"),intern = T)
	file.remove(filesIn[i,"ft"])
	filesIn[i,]
}

print(paste(date(),"Trimming and QC checking"))
cat("",file=file.path(ED,"exclude"))
filesInEnvironment<-foreach(i=which(filesIn[,"gr"]=="environment"),.combine = rbind) %dopar% filesInGen(i)
RevCompl <- function(x="") chartr("ATGCatgcUu","TACGtacgAa",sapply(lapply(strsplit(x,NULL),rev),paste,collapse=""))
exclude<-c()
for(f in filesInEnvironment[filesInEnvironment[,"gr"]=="environment","name"]){
	tmp<-read.table(paste0(file.path(ED,f,f),".reads.gz"))
	exclude<-c(exclude,tmp[tmp[,1]>1,2])
}
if(length(exclude)>0){
	exclude<-unique(exclude)
	exclude<-c(exclude,unlist(lapply(exclude, RevCompl)))
	exclude<-exclude[order(exclude)]
	cat(unique(exclude),file=file.path(ED,"exclude"), sep="\n", append=F)
}

filesIn<-foreach(i=which(filesIn[,"gr"]!="environment"),.combine = rbind) %dopar% filesInGen(i)
filesIn<-rbind(filesIn,filesInEnvironment)
zip(file.path(ED,paste0(Exp,"_fastQC.zip")),files=dir(file.path(ED,"qc"),".html",full.names = T),flags="-joq9")

setwd(file.path(ED,"qc"))
con<-file(file.path("..",paste0(Exp,"_fastQC.html")),"wt")
index<-dir(".",".html$")
cat(paste0('<html><head><link rel="stylesheet" type="text/css" href="/shared/shiny.css"/></head><body><table width="100%" border="2">',
		'<tr><td width="10%">File</td><td rowspan="',length(index),' width="90%""><iframe srcdoc="" name="QC" width="100%" height="',length(index)*40,'px"> </iframe></td></tr>'),file = con)
for(i in index) cat(paste0("<tr align=center><td><a href=qc/",i," target='QC'>",i,"</a></td></tr>"),file = con)
cat(paste("</table></body></html>"),file = con)
close(con)
setwd(wd)

save(filesIn,file = file.path(ED,"filesIn.RData"))
