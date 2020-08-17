library(shiny)
library(gdata)
library(foreach)
library(doMC)
core <- parallel::detectCores(all.tests =TRUE,logical =FALSE)-1
registerDoMC(cores = core)

# library(promises)
# library(future)
# plan("multiprocess")

wd<-sub("/shiny$","",getwd())
setwd(wd)
source(file.path(wd,"bin","utils.R"))

# examplesdir<-file.path(wd,"www","upload","example-samples")
# if(!dir.exists(examplesdir)) dir.create(examplesdir,recursive = TRUE)
# examples<-dir(path=examplesdir,pattern = ".(fastq|fq|fasta|fa|bam|sam|cram)(|.(gz|bz2|xz))$", full.names = FALSE, recursive = TRUE, include.dirs = TRUE)
# if(length(examples)>0){
#   examples<-cbind(file=file.path("example-samples",examples),
# 				size=humanReadable(file.info(file.path(examplesdir,examples))$size),
# 				date=format(file.info(file.path(examplesdir,examples))$mtime,"%d.%m.%Y %H:%M:%OS"))
# } else examples<-rbind(rep(NA,3))[-1,]
# if(!dir.exists(file.path(wd,"www","upload"))) dir.create(file.path(wd,"www","upload"),recursive = TRUE)
# serverFiles<-dir(path = file.path("www","upload"),pattern = ".(fastq|fq|fasta|fa|bam|sam|cram)(|.(gz|bz2|xz))$",full.names = FALSE, recursive = TRUE, include.dirs = TRUE)
# if(length(grep("^example-samples",serverFiles))>0) serverFiles<-serverFiles[-grep("^example-samples",serverFiles)]
# if(length(serverFiles)>0){
# 	serverFiles<-cbind(file=serverFiles,
# 				size=humanReadable(file.info(file.path("www","upload",serverFiles))$size),
# 				date=format(file.info(file.path("www","upload",serverFiles))$mtime,"%d.%m.%Y %H:%M:%OS"))
# } else serverFiles<-rbind(rep(NA,3))[-1,]
filesUploaded           <- rbind(rep(NA,3))[-1,]
colnames(filesUploaded) <- c("file","size","date") #  <- colnames(examples) <- colnames(serverFiles)
# groups           <- rbind(rep(NA,6))[-1,]
# colnames(groups) <- c("file","size","date","test","control","ignore")
if(!dir.exists(file.path(wd,"www","db"))) dir.create(file.path(wd,"www","db"),recursive = TRUE)
if(!exists("FilesIn"))
	if(file.exists(file.path(wd,"www","db","FilesIn.RData"))){
		load(file.path(wd,"www","db","FilesIn.RData"))
	} else { FilesIn <- rbind(rep(NA,3))[-1,]; colnames(FilesIn) <- c("file","size","date"); }
if(!exists("GroupsSel"))
	if(file.exists(file.path(wd,"www","db","GroupsSel.RData"))){
		load(file.path(wd,"www","db","GroupsSel.RData"))
	} else { GroupsSel <- rbind(rep(NA,6))[-1,]; colnames(GroupsSel) <- c("file","size","date","test","control","ignore"); }

species<-c("homo_sapiens","mus_musculus")
if(!exists("tsize")){
	specie<-"homo_sapiens"
	tsize<-"2000"
	Rep<-2
	blast<-"nr"
	qc<-20
	ad3<-"TGGAATTCTCGGGTGCCAAGG # Illumina TruSeq Small RNA"  # ad3<-"ATCACCGACTGCCCATAGAGAG"  # Ion Torrent
	ad5<-"GTTCAGAGTTCTACAGTCCGACGATC # Illumina TruSeq Small RNA" # ad5<-"CCAAGGCG"  # Ion Torrent
	sizerange<-c(10,42)
	lim<-2
	log2FoldChange<-1
	padj<-0.05
	email<-""
	smtpServer<-""
	if(file.exists(file.path(wd,"www","db","Config.RData"))){
		load(file.path(wd,"www","db","Config.RData"))
	} else {
		Exp <- ""
	}
}


