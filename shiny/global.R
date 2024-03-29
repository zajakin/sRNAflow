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
Sys.umask("0")
Sys.setlocale(category = "LC_ALL", locale = "en_US.utf8")

filesUploaded           <- rbind(rep(NA,3))[-1,]
colnames(filesUploaded) <- c("file","size","date")
if(!dir.exists(file.path(wd,"www","db"))){
	dir.create(file.path(wd,"www","db"),recursive = TRUE, mode = "0777")
	examplesdir<-file.path(wd,"www","results","Example")
	if(!dir.exists(examplesdir)) dir.create(examplesdir,recursive = TRUE, mode = "0777")
	file.copy(dir(file.path(wd,"examples","results"),full.names =TRUE),examplesdir, recursive = TRUE, copy.mode = FALSE, copy.date = TRUE)
	unzip(file.path(examplesdir,"Example_sRNAflow_diagrams.zip"),exdir = file.path(examplesdir,"species_diagrams"))
	unzip(file.path(examplesdir,"Example_fastQC.zip"),exdir = file.path(examplesdir,"qc"))
	examplesdir<-file.path(wd,"www","upload","example-samples")
	if(!dir.exists(examplesdir)) dir.create(examplesdir,recursive = TRUE, mode = "0777")
	file.copy(dir(file.path(wd,"examples","samples"),full.names =TRUE),examplesdir, recursive = FALSE, copy.mode = FALSE, copy.date = TRUE)
}
if(!exists("FilesIn"))
	if(file.exists(file.path(wd,"www","db","FilesIn.RData"))){
		load(file.path(wd,"www","db","FilesIn.RData"))
	} else { FilesIn <- rbind(rep(NA,3))[-1,]; colnames(FilesIn) <- c("file","size","date"); }
if(!exists("GroupsSel"))
	if(file.exists(file.path(wd,"www","db","GroupsSel.RData"))){
		load(file.path(wd,"www","db","GroupsSel.RData"))
	} else { GroupsSel <- rbind(rep(NA,7))[-1,]; colnames(GroupsSel) <- c("file","size","date","test","control","environment","ignore"); }

species<-c("homo_sapiens")  #,"mus_musculus")
if(!exists("tsize")){
	specie<-"homo_sapiens"
	strategy<-"metagenome"
	tsize<-"200"
	Rep<-1
	blast<-"main specie & bacteria+"
	qc<-8
	ad3<-"TGGAATTCTCGGGTGCCAAGG # Illumina TruSeq Small RNA"  # ad3<-"ATCACCGACTGCCCATAGAGAG"  # Ion Torrent # ad3<-"AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCGGCTTG" # universal Illumina adapter
	ad5<-"GTTCAGAGTTCTACAGTCCGACGATC # Illumina TruSeq Small RNA" # ad5<-"CCAAGGCG"  # Ion Torrent
	sizerange<-c(18,50)
	lim<-2
	limS<-1
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


