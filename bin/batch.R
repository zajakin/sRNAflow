#!/usr/bin/R --no-save
wd<-sub("/bin$","",getwd())
setwd(wd)
source("shiny/global.R")
#TODO  GUI for FilesIn; GroupsSel; Exp; Rep

# apt install ncbi-blast+ ncbi-entrez-direct # blastn –taxidlist meta.txids
system(paste("awk -F'\t' '{print $2}' bin/taxids_for_blast.tsv | xargs -l /usr/bin/get_species_taxids -t > db/meta.txids"))
specie<-"homo_sapiens"
tsize<-"2000"
Rep<-2
lim<-2
log2FoldChange<-1
padj<-0.05

save(FilesIn,GroupsSel,Exp,Rep,specie,tsize,Rep,lim,log2FoldChange,padj,file = file.path(wd,"data",Exp,"settings.RData"))
load(file.path(wd,"data",Exp,"settings.RData"))

source("bin/sRNAflow_filesIn_subsets.R")
load(file.path(wd,"data",Exp,"filesIn.RData"))

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
load("db/genomes/genomesdb.RData")
#download genomes  index generation   #####
source("bin/sRNAflow_downloadGenomes.R")

system(paste("ln -fs ",file.path(wd,"db","gtf_biotypes"),file.path(wd,"data",Exp,"genomes","gtf_biotypes")))
system(paste("ln -fs ",file.path(wd,"db","genomes","homo_sapiens.gtf"),file.path(wd,"data",Exp,"genomes","genomes.gtf")))

#mapping
library(foreach)
library(doMC)
registerDoMC()
err<-foreach(file=filesIn[,"wf"]) %dopar% {
	system(paste0(wd,"/bin/sRNAflow_mapping_and_RNA_catalog.sh ",file.path(wd,file)," ",file.path(wd,"data",Exp)),intern = TRUE)
}
#United table
system(paste0(wd,"/bin/sRNAflow_united_table_of_mapping_and_RNA_catalogs.sh ",file.path(wd,"data",Exp)),intern = TRUE)

source("bin/sRNAflow_DESeq2.R")

#TODO isomiR-SEA (?)

#TODO miRanda,TarBase, MiRtaget2 (?)

#TODO MultiQC (?)

#TODO mail

