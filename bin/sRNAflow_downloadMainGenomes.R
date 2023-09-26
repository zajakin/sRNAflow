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
	}
}
print(paste(date(),"Download main genome files"))
timeout<-getOption('timeout')
options(timeout=999999)
# as.vector(md5sum(dir(R.home(), pattern = "^COPY", full.names = TRUE)))

#### Pull taxonomy from BLAST
if(!file.exists(file.path(wd,"www","db","blast.txids")) || difftime(Sys.time(),file.mtime(file.path(wd,"www","db","blast.txids")),units = "days")>30 || 
   difftime(file.mtime(file.path(wd,"bin","taxids_for_blast.tsv")),file.mtime(file.path(wd,"www","db","blast.txids")),units="secs")>0)
	system(paste("gawk -F'\t' '{print $2}'",file.path(wd,"bin","taxids_for_blast.tsv"),"| xargs --max-procs=1 -l get_species_taxids -t >",
				 file.path(wd,"www","db","blast.txids")))
blasttaxids<-read.table(file.path(wd,"www","db","blast.txids"))[,1]

#### Pull taxonomy from Krona  # system("kronatools_updateTaxonomy")
if(!dir.exists(file.path(wd,"www","db","taxonomy"))) dir.create(file.path(wd,"www","db","taxonomy"),recursive = TRUE, mode = "0777")
system(paste(file.path("","usr","bin","kronatools_updateTaxonomy && cp"),file.path("","var","lib","kronatools","taxonomy","*"),file.path(wd,"www","db","taxonomy")))

#### Make meta.txids from taxonomy.tab
if(!file.exists(file.path(wd,"www","db","meta.txids")) 
   || difftime(file.mtime(file.path(wd,"bin","taxids_for_blast.tsv")),file.mtime(file.path(wd,"www","db","meta.txids")),units="secs")>0
   || difftime(file.mtime(file.path(wd,"www","db","blast.txids")),file.mtime(file.path(wd,"www","db","meta.txids")),units="secs")>0
   || difftime(file.mtime(file.path(wd,"www","db","taxonomy","taxonomy.tab")),file.mtime(file.path(wd,"www","db","meta.txids")),units="secs")>0){
		taxids<-read.table(file.path(wd,"bin","taxids_for_blast.tsv"),sep = "\t")
		taxtab<-read.table(file.path(wd,"www","db","taxonomy","taxonomy.tab"),sep = "\t",quote = "")
		taxout<-as.character(unique(taxids[,2]))
		j<-as.character(taxtab[as.character(taxtab[,3]) %in% taxout,1])
		while(length(j)>0){
			taxout<-c(taxout,j)
			j<-as.character(taxtab[as.character(taxtab[,3]) %in% j,1])
		}
		taxout<-unique(c(taxout,blasttaxids))
		cat(taxout[order(taxout)],file=file.path(wd,"www","db","meta.txids"),sep="\n")
}

if(!dir.exists(file.path(wd,"www","db","genomes"))) dir.create(file.path(wd,"www","db","genomes"),recursive = TRUE, mode = "0777")

i<-file.path(wd,"www","db","genomes",paste0(specie,".gtf"))
if(!file.exists(paste0(i,".gz")) || difftime(Sys.time(),file.mtime(paste0(i,".gz")),units = "days")>30){
	getDBfile(c('ftp://ftp.ensembl.org/pub/current_gtf/',specie),specie,".gtf.gz",".gtf.gz")
	tax<-system(paste0("gawk -F'\t' 'tolower($5) ~/^",sub("_"," ",specie),"$/{print $1}' www/db/taxonomy/taxonomy.tab"),intern = TRUE)
	system(paste0("cat ",i," | sed -E '/(^#|^$)/!s/^/",tax,"_",specie,"_/' > ",sub(".gtf","_tax.gtf",i)),intern = TRUE)
	system(paste("pigz -9f",i,sub(".gtf","_tax.gtf",i)))
}

if(!file.exists(file.path(wd,"www","db","genomes",paste0(specie,".fa"))) || difftime(Sys.time(),file.mtime(file.path(wd,"www","db","genomes",paste0(specie,".fa"))),units = "days")>30){
	getDBfile(c('ftp://ftp.ensembl.org/pub/current_fasta/',specie,'/dna'),specie,".dna.primary_assembly.fa.gz",".fa.gz")
	dir.create(file.path(wd,"www","db","genomes","bowtie2",specie),recursive = TRUE, mode = "0777")
	system(paste0("bowtie2-build --threads ",core," ",file.path(wd,"www","db","genomes",specie),".fa ",file.path(wd,"www","db","genomes","bowtie2",specie,specie)," > ",file.path(wd,"www","db","genomes","bowtie2",specie,"bowtie2-build.log")))
# system(paste0("zcat ",specie,".dna.primary_assembly.fa.gz > ",specie,".fa"))
# system(paste0("zcat ",specie,".gtf.gz > ",specie,".gtf"))
}
if(!file.exists(file.path(wd,"www","db","genomes","univec.fa")) || difftime(Sys.time(),file.mtime(file.path(wd,"www","db","genomes","univec.fa")),units = "days")>30){
	download.file("https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core",file.path("www","db","genomes","univec.fa"))
	dir.create(file.path(wd,"www","db","genomes","bowtie2","univec"),recursive = TRUE, mode = "0777")
	system(paste("bowtie2-build --threads",core,file.path(wd,"www","db","genomes","univec.fa"),file.path(wd,"www","db","genomes","bowtie2","univec","univec"),">",file.path(wd,"www","db","genomes","bowtie2","univec","bowtie2-build.log")))
}
if(!file.exists(file.path(wd,"www","db","genomes","mature.fa")) || difftime(Sys.time(),file.mtime(file.path(wd,"www","db","genomes","mature.fa")),units = "days")>30){
	download.file("https://mirbase.org/download/mature.fa",file.path("www","db","genomes","mature.fa"))
#	system(paste("pigz -df",file.path("www","db","genomes","mature.fa.gz")))
# pigz -cd $DV.fa.gz | fasta_formatter | sed '/^[^>]/ y/uU/tT/' > $DV.fa
}
if(!file.exists(file.path(wd,"www","db","genomes","ensemblgenomes.txt")) || difftime(Sys.time(),file.mtime(file.path(wd,"www","db","genomes","ensemblgenomes.txt")),units = "days")>30){
	download.file("ftp://ftp.ensemblgenomes.org/pub/current/species.txt",file.path("www","db","genomes","ensemblgenomes.txt"))
	system("sed -i '1 s/$/\tNA/' www/db/genomes/ensemblgenomes.txt",intern = TRUE)
}
if(!file.exists(file.path(wd,"www","db","genomes","genbank.txt")) || difftime(Sys.time(),file.mtime(file.path(wd,"www","db","genomes","genbank.txt")),units = "days")>30)
	download.file("https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt","www/db/genomes/genbank.txt")
if(!file.exists(file.path(wd,"www","db","genomes","refseq.txt")) || difftime(Sys.time(),file.mtime(file.path(wd,"www","db","genomes","refseq.txt")),units = "days")>30)
	download.file("https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt","www/db/genomes/refseq.txt")

if(!dir.exists(file.path(wd,"www","db","gtf_biotypes"))){
	dir.create(file.path(wd,"www","db","gtf_biotypes"),recursive = TRUE, mode = "0777")
	file.copy(dir(file.path(wd,"gtf_biotypes"),full.names =TRUE),file.path(wd,"www","db","gtf_biotypes"), recursive = FALSE, copy.mode = FALSE, copy.date = TRUE)
}

options(timeout=timeout)

