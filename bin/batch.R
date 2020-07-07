#!/usr/bin/R --no-save
wd<-sub("/bin$","",getwd())
setwd(wd)
source("shiny/global.R")
#TODO  GUI for FilesIn; GroupsSel; Exp; Rep
specie<-"homo_sapiens"
tsize<-"2000"
Rep<-2
lim<-2
log2FoldChange<-1
padj<-0.05

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
		filesIn[i,paste0("ft",r)]<-paste0(d,"faTab/",s,"_random",tsize,".",r,".faTab")
	}
}

save(filesIn,FilesIn,GroupsSel,Exp,Rep,specie,tsize,Rep,lim,log2FoldChange,padj,file = file.path(wd,"data",Exp,"settings.RData"))
load(file.path(wd,"data",Exp,"settings.RData"))

#blast per sample
#args <- c("R","--no-save",paste0("data/test"),paste0("../example-samples/Cancer_1.fa"),4)
for(idr in 1:nrow(filesIn)){
	for(re in 1:as.numeric(Rep)){
		args <- c("sRNAflow","--no-save",file.path(wd,filesIn[idr,"wd"]),paste0(filesIn[idr,"name"],"_random",tsize,".",re),file.path(wd,filesIn[idr,paste0("ft",re)]),4)
		source(file.path(wd,"bin/sRNAflow_blast_per_sample.R"))
	}
}

#aggregate species  ########
source("bin/sRNAflow_aggregate_species.R")

#download main genomes from Ensembl  ########
source("bin/sRNAflow_downloadMainGenomes.R")

#create GTF files ####


#download genomes  index generation   #####
source("bin/sRNAflow_downloadGenomes.R")

#TODO download --- genomes.gtf ####


# system(paste("sed -i -E '/(^#|^$)/!s/^/9606_homo_sapiens_/' db/genomes/homo_sapiens.gtf"),intern = TRUE)
# system(paste("sed -i -E '/(^#|^$)/!s/^/9606_homo_sapiens_/' db/gtf_biotypes/*.gtf"),intern = TRUE)
system(paste("ln -fs ",file.path(wd,"db","gtf_biotypes"),file.path(wd,"data",Exp,"genomes","gtf_biotypes")))
system(paste("ln -fs ",file.path(wd,"db","genomes","homo_sapiens.gtf"),file.path(wd,"data",Exp,"genomes","genomes.gtf")))

#mapping
for(file in filesIn[,"wf"]) system(paste0(wd,"/bin/sRNAflow_mapping_and_RNA_catalog.sh ",file.path(wd,file)," ",file.path(wd,"data",Exp)),intern = TRUE)

#United table
system(paste0(wd,"/bin/sRNAflow_united_table_of_mapping_and_RNA_catalogs.sh ",file.path(wd,"data",Exp)),intern = TRUE)

source("bin/sRNAflow_DESeq2.R")

#TODO isomiR-SEA (?)

#TODO miRanda,TarBase, MiRtaget2 (?)

#TODO MultiQC (?)

#TODO mail

