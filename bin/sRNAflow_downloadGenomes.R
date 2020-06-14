#!/usr/bin/env Rscript --vanilla 
WD<-file.path("data",Exp,"genomes")
if(!dir.exists(file.path(WD,"archive"))) dir.create(file.path(WD,"archive"),recursive = TRUE)
i<-0
ext<-"_genomic.fna.gz"
fa<-file.path(WD,"genomes.fa")
file.create(fa)
for(id in species99[,"id"]){
	print(paste(i<-i+1,id,species99[as.character(id),2],date()))
	spname<-sub(" .*","",sub(" ","_",species99[as.character(id),2]))
	if(id=="9606"){
		system(paste0("cat db/genomes/",specie,".fa | sed 's/^>/>",id,"_",specie,"_/g' >> ",fa))
		next
	} 
	urls<-refseq[refseq[,"taxid"]==id,"ftp_path"]
	if(length(urls)==0) urls<-refseq[refseq[,"species_taxid"]==id,"ftp_path"]
	if(length(urls)==0) urls<-genbank[genbank[,"taxid"]==id,"ftp_path"]
	if(length(urls)==0) urls<-genbank[genbank[,"species_taxid"]==id,"ftp_path"]
	if(length(urls)==0) urls<-ensemblgenomes[ensemblgenomes[,"taxonomy_id"]==id,"ftp_path"]
	dir.create(file.path(WD,"archive",id))
	for(url in urls[1]){
		if(!file.exists(file.path(WD,"archive",id,paste0(basename(url),ext))))
			download.file(paste0(url,"/",basename(url),ext),file.path(WD,"archive",id,paste0(basename(url),ext)))
		system(paste0("gunzip -c ",file.path(WD,"archive",id,paste0(basename(url),ext))," | sed 's/^>/>",id,"_",spname,"_/g' >> ",fa)  )
	}
}

if(file.size(fa)>(2^32-1)){ system(paste("bowtie2-build --threads ",core," --large-index ",fa,file.path(WD,"genomes")))
} else system(paste("bowtie2-build --threads ",core,fa,file.path(WD,"genomes")))

# as.vector(md5sum(dir(R.home(), pattern = "^COPY", full.names = TRUE)))
