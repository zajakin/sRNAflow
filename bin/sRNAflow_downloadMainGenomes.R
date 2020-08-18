#!/usr/bin/env Rscript --vanilla 
getDir<-function(base,sep=''){
	filenames<-c()
	if(zz<- file(paste(c(base,'/'),collapse=""),open="r",method="libcurl")){
		if(substr(base[1],1,5)=='http:') sep='"'
		filenames<-scan(zz,what=character(),sep=sep)
		close(zz)
	}
	return(filenames)
}
getDBfile<-function(base=c('http://ftp.ensembl.org/pub/current_gtf/',specie),sp=specie,ext=".gtf.gz",sep='',path=file.path("www","db","genomes")){
	if(length(filenames<-getDir(base))>0){
		filenames <- as.character(filenames[grep(ext,filenames)])
		tmp<-sapply(filenames,nchar)
		filenames <-filenames[tmp==min(tmp)][1]
		download.file(paste(c(base,'/',filenames),collapse=""),file.path(path,paste0(sp,ext)),"auto",mode = "wb")
		ext2<-gsub(".*\\.",".",sub(".gz","",ext))
		out<-file(file.path(path,paste0(sp,ext2)),"wt")
		zz <-gzcon(file(file.path(path,paste0(sp,ext)),"r"))
		tmp<-readLines(zz)
		close(zz)
		writeLines(tmp,out,sep="\n")
		close(out)
		file.remove(file.path(path,paste0(sp,ext)))
	}
}

getDBfile(c('ftp://ftp.ensembl.org/pub/current_gtf/',specie),specie,".gtf.gz")
getDBfile(c('ftp://ftp.ensembl.org/pub/current_fasta/',specie,'/dna'),specie,".dna.primary_assembly.fa.gz")
# system(paste0("zcat ",specie,".dna.primary_assembly.fa.gz > ",specie,".fa"))
# system(paste0("zcat ",specie,".gtf.gz > ",specie,".gtf"))

download.file("ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz",file.path("www","db","genomes","mature.fa.gz"))
system(paste("pigz -d",file.path("www","db","genomes","mature.fa.gz")))
# pigz -cd $DV.fa.gz | fasta_formatter | sed '/^[^>]/ y/uU/tT/' > $DV.fa

download.file("ftp://ftp.ensemblgenomes.org/pub/current/species.txt",file.path("www","db","genomes","ensemblgenomes.txt"))
system("sed -i '1 s/$/\tNA/' www/db/genomes/ensemblgenomes.txt",intern = TRUE)
ensemblgenomes<-read.table("www/db/genomes/ensemblgenomes.txt", comment.char="",skip=0,quote = "", header = TRUE, sep = "\t",dec = ".", na.strings = "",as.is = FALSE)
head(ensemblgenomes) # ,row.names = NULL 
dim(ensemblgenomes)

download.file("ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt","www/db/genomes/genbank.txt")
genbank<-read.table("www/db/genomes/genbank.txt", comment.char="",skip=1,quote = "", header = TRUE, sep = "\t",dec = ".", na.strings = "",as.is = TRUE)
# rownames(genbank)<-genbank[,"taxid"]
head(genbank)
dim(genbank)

download.file("ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt","www/db/genomes/refseq.txt")
refseq<-read.table("www/db/genomes/refseq.txt", comment.char="",skip=1,quote = "", header = TRUE, sep = "\t",dec = ".", na.strings = "",as.is = TRUE)
# rownames(refseq)<-refseq[,"taxid"]
head(refseq)
dim(refseq)

save(ensemblgenomes,genbank,refseq,file="www/db/genomes/genomesdb.RData")

#NOTE create GTF files, remove united ####
#NOTE download --- genomes.gtf ####
system(paste("sed -i -E '/(^#|^$)/!s/^/9606_homo_sapiens_/' www/db/genomes/homo_sapiens.gtf"),intern = TRUE)
# system(paste("sed -i -E '/(^#|^$)/!s/^/9606_homo_sapiens_/' www/db/gtf_biotypes/*.gtf"),intern = TRUE)

system(paste("ln -fs ",file.path(wd,"www","db","gtf_biotypes"),file.path(ED,"genomes","gtf_biotypes")))
system(paste("ln -fs ",file.path(wd,"www","db","genomes","homo_sapiens.gtf"),file.path(ED,"genomes","genomes.gtf")))

