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
Sys.setlocale(category = "LC_ALL", locale = "en_US.utf8")

filesUploaded           <- rbind(rep(NA,3))[-1,]
colnames(filesUploaded) <- c("file","size","date")
if(!dir.exists(file.path(wd,"www","db"))) dir.create(file.path(wd,"www","db"),recursive = TRUE, mode = "0777")
if(!exists("FilesIn"))
	if(file.exists(file.path(wd,"www","db","FilesIn.RData"))){
		load(file.path(wd,"www","db","FilesIn.RData"))
	} else { FilesIn <- rbind(rep(NA,3))[-1,]; colnames(FilesIn) <- c("file","size","date"); }
if(!exists("GroupsSel"))
	if(file.exists(file.path(wd,"www","db","GroupsSel.RData"))){
		load(file.path(wd,"www","db","GroupsSel.RData"))
	} else { GroupsSel <- rbind(rep(NA,6))[-1,]; colnames(GroupsSel) <- c("file","size","date","test","control","ignore"); }

species<-c("homo_sapiens")  #,"mus_musculus")
if(!exists("tsize")){
	specie<-"homo_sapiens"
	strategy<-"metagenome"
	tsize<-"200"
	Rep<-2
	blast<-"nr/nt"
	qc<-20
	ad3<-"TGGAATTCTCGGGTGCCAAGG # Illumina TruSeq Small RNA"  # ad3<-"ATCACCGACTGCCCATAGAGAG"  # Ion Torrent # ad3<-"AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCGGCTTG" # universal Illumina adapter
	ad5<-"GTTCAGAGTTCTACAGTCCGACGATC # Illumina TruSeq Small RNA" # ad5<-"CCAAGGCG"  # Ion Torrent
	sizerange<-c(18,50)
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


