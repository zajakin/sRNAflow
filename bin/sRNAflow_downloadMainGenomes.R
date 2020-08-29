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
getDBfile<-function(base=c('http://ftp.ensembl.org/pub/current_gtf/',specie),sp=specie,ext=".gtf.gz",ext2=".gtf.gz",sep='',path=file.path("www","db","genomes")){
	if(length(filenames<-getDir(base))>0){
		filenames <- as.character(filenames[grep(ext,filenames)])
		tmp<-sapply(filenames,nchar)
		filenames <-filenames[tmp==min(tmp)][1]
		download.file(paste(c(base,'/',filenames),collapse=""),file.path(path,paste0(sp,ext2)),"auto",mode = "wb")
		system(paste("pigz -df",file.path(path,paste0(sp,ext2))))
		# out<-file(file.path(path,paste0(sp,ext2)),"wt")
		# zz <-gzcon(file(file.path(path,paste0(sp,ext)),"r"))
		# tmp<-readLines(zz)
		# close(zz)
		# writeLines(tmp,out,sep="\n")
		# close(out)
		# file.remove(file.path(path,paste0(sp,ext)))
	}
}

# as.vector(md5sum(dir(R.home(), pattern = "^COPY", full.names = TRUE)))
if(!dir.exists(file.path(wd,"www","db","taxonomy"))) dir.create(file.path(wd,"www","db","taxonomy"),recursive = TRUE, mode = "0777")
system(paste(file.path(wd,"Krona","KronaTools","updateTaxonomy.sh"),file.path(wd,"www","db","taxonomy")))
# system("kronatools_updateTaxonomy")

if(!dir.exists(file.path(wd,"www","db","genomes"))) dir.create(file.path(wd,"www","db","genomes"),recursive = TRUE, mode = "0777")

if(!file.exists(file.path(wd,"www","db","genomes",paste0(specie,".gtf")))){
	getDBfile(c('ftp://ftp.ensembl.org/pub/current_gtf/',specie),specie,".gtf.gz",".gtf.gz")
	tax<-system(paste0("gawk -F'\t' 'tolower($5) ~/^",sub("_"," ",specie),"$/{print $1}' www/db/taxonomy/taxonomy.tab"),intern = TRUE)
	system(paste0("sed -i -E '/(^#|^$)/!s/^/",tax,"_",specie,"_/' www/db/genomes/",specie,".gtf"),intern = TRUE)
}
if(!file.exists(file.path(wd,"www","db","genomes",paste0(specie,".fa")))) getDBfile(c('ftp://ftp.ensembl.org/pub/current_fasta/',specie,'/dna'),specie,".dna.primary_assembly.fa.gz",".fa.gz")
# system(paste0("zcat ",specie,".dna.primary_assembly.fa.gz > ",specie,".fa"))
# system(paste0("zcat ",specie,".gtf.gz > ",specie,".gtf"))

if(!file.exists(file.path(wd,"www","db","genomes","mature.fa"))){
	download.file("ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz",file.path("www","db","genomes","mature.fa.gz"))
	system(paste("pigz -df",file.path("www","db","genomes","mature.fa.gz")))
# pigz -cd $DV.fa.gz | fasta_formatter | sed '/^[^>]/ y/uU/tT/' > $DV.fa
}
if(!file.exists(file.path(wd,"www","db","genomes","ensemblgenomes.txt"))){
	download.file("ftp://ftp.ensemblgenomes.org/pub/current/species.txt",file.path("www","db","genomes","ensemblgenomes.txt"))
	system("sed -i '1 s/$/\tNA/' www/db/genomes/ensemblgenomes.txt",intern = TRUE)
}
if(!file.exists(file.path(wd,"www","db","genomes","genbank.txt"))) download.file("ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt","www/db/genomes/genbank.txt")
if(!file.exists(file.path(wd,"www","db","genomes","refseq.txt"))) download.file("ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt","www/db/genomes/refseq.txt")

if(!file.exists(file.path(wd,"www","db","meta.txids")))
	system(paste("gawk -F'\t' '{print $2}'",file.path(wd,"bin","taxids_for_blast.tsv"),"| xargs -l get_species_taxids -t >",file.path(wd,"www","db","meta.txids")))

if(!dir.exists(file.path(wd,"www","db","gtf_biotypes"))){
	dir.create(file.path(wd,"www","db","gtf_biotypes"),recursive = TRUE, mode = "0777")
	file.copy(dir(file.path(wd,"gtf_biotypes"),full.names =TRUE),file.path(wd,"www","db","gtf_biotypes"), recursive = FALSE, copy.mode = FALSE, copy.date = TRUE)
}

