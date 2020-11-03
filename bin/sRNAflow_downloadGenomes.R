#!/usr/bin/env Rscript --vanilla 
WD<-file.path(ED,"genomes")
archive<-file.path(wd,"www","db","genomes","archive")
ensemblgenomes<-read.table("www/db/genomes/ensemblgenomes.txt", comment.char="",skip=0,quote = "", header = TRUE, sep = "\t",dec = ".", na.strings = "",as.is = FALSE)
genbank<-read.table("www/db/genomes/genbank.txt", comment.char="",skip=1,quote = "", header = TRUE, sep = "\t",dec = ".", na.strings = "",as.is = TRUE)
refseq<-read.table("www/db/genomes/refseq.txt", comment.char="",skip=1,quote = "", header = TRUE, sep = "\t",dec = ".", na.strings = "",as.is = TRUE)
# rownames(genbank)<-genbank[,"taxid"]
# rownames(refseq)<-refseq[,"taxid"]
# save(ensemblgenomes,genbank,refseq,file="www/db/genomes/genomesdb.RData")
# load(file.path("www","db","genomes","genomesdb.RData"))
if(!dir.exists(WD)) dir.create(WD,recursive = TRUE)
if(!dir.exists(archive)) dir.create(archive,recursive = TRUE)
i<-0
ext<-"_genomic.fna.gz"
fa<-file.path(WD,"genomes.fa")
file.create(fa)
timeout<-getOption('timeout')
options(timeout=999999)
for(id in species99[,"id"]){
	print(paste(i<-i+1,id,species99[as.character(id),2],date()))
	spname<-sub(" .*","",sub(" ","_",species99[as.character(id),2]))
	if(id=="9606"){
		system(paste0("cat www/db/genomes/",specie,".fa | sed 's/^>/>",id,"_",specie,"_/g' >> ",fa))
		next
	} 
	urls<-refseq[refseq[,"taxid"]==id,"ftp_path"]
	if(length(urls)==0) urls<-refseq[refseq[,"species_taxid"]==id,"ftp_path"]
	if(length(urls)==0) urls<-genbank[genbank[,"taxid"]==id,"ftp_path"]
	if(length(urls)==0) urls<-genbank[genbank[,"species_taxid"]==id,"ftp_path"]
	if(length(urls)==0) urls<-ensemblgenomes[ensemblgenomes[,"taxonomy_id"]==id,"ftp_path"]
	dir.create(file.path(archive,id))
	for(url in urls[1]){
		if(!file.exists(file.path(archive,id,paste0(basename(url),ext)))){
			download.file(paste0(url,"/",basename(url),ext),file.path(archive,id,paste0(basename(url),ext)))
		}
		system(paste0("pigz -cd ",file.path(archive,id,paste0(basename(url),ext))," | sed 's/^>/>",id,"_",spname,"_/g' >> ",fa)  )
	}
	#IDEA Select representative contig per specie (by clustering?)
}
options(timeout=timeout)
rm(ensemblgenomes,genbank,refseq)

if(file.size(fa)>(2^32-1)){ system(paste("bowtie2-build --threads ",core," --large-index ",fa,file.path(WD,"genomes"),">",file.path(WD,"bowtie2-build.log")))
} else system(paste("bowtie2-build --threads ",core,fa,file.path(WD,"genomes"),">",file.path(WD,"bowtie2-build.log")))

system(paste("ln -fs ",file.path(wd,"www","db","gtf_biotypes"),WD))
system(paste("ln -fs ",file.path(wd,"www","db","genomes","homo_sapiens.gtf"),file.path(WD,"genomes.gtf")))
