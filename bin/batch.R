#!/usr/bin/R --no-save
# args <- commandArgs(TRUE)
wd<-sub("/bin$","",getwd())
setwd(wd)
source("shiny/global.R")
# apt install ncbi-blast+ ncbi-entrez-direct # blastn –taxidlist meta.txids
# system(paste("awk -F'\t' '{print $2}' bin/taxids_for_blast.tsv | xargs -l /usr/bin/get_species_taxids -t > www/db/meta.txids"))
ED<-file.path(wd,"www",Exp)
dir.create(ED,recursive = TRUE)
save(FilesIn,GroupsSel,Exp,specie,tsize,Rep,blast,ad3,ad5,sizerange,lim,log2FoldChange,padj,email,smtpServer,file = file.path(ED,"settings.RData"))
load(file.path(ED,"settings.RData"))

source("bin/sRNAflow_filesIn_subsets.R")
load(file.path(ED,"filesIn.RData"))

source(file.path(wd,"bin","sRNAflow_blast_per_sample.R"))
err<-foreach(idr=1:nrow(filesIn),.verbose = T) %dopar% for(re in 1:as.numeric(Rep)) blast_per_sample(idr,re,wd,filesIn,tsize,core)

#aggregate species  ########
source("bin/sRNAflow_aggregate_species.R")

#download main genomes from Ensembl  ########
source("bin/sRNAflow_downloadMainGenomes.R")
load(file.path("www","db","genomes","genomesdb.RData"))
#download genomes  index generation   #####
source("bin/sRNAflow_downloadGenomes.R")

system(paste("ln -fs ",file.path(wd,"www","db","gtf_biotypes"),file.path(ED,"genomes","gtf_biotypes")))
system(paste("ln -fs ",file.path(wd,"www","db","genomes","homo_sapiens.gtf"),file.path(ED,"genomes","genomes.gtf")))

#mapping o:f:r:l:a
err<-foreach(i=1:nrow(filesIn)) %dopar% {
	system(paste(file.path(wd,"bin","sRNAflow_mapping_and_RNA_catalog.sh"),
				 "-n",filesIn[i,"name"],
				 "-r",file.path(wd,filesIn[i,"rf"]),
				 "-f",file.path(wd,filesIn[i,"wf"]),
				 "-t",filesIn[i,"type"],
				 "-o",ED),intern = TRUE)
}
#United table
system(paste0(wd,"/bin/sRNAflow_united_table_of_mapping_and_RNA_catalogs.sh ",ED),intern = TRUE)

source("bin/sRNAflow_DESeq2.R")

# OPTIMIR https://github.com/FlorianThibord/OptimiR

# https://www.bioconductor.org/packages/release/bioc/html/isomiRs.html
# BiocManager::install("isomiRs") 
# apt install mirtop libpng-dev
# library("isomiRs")
# data(mirData)
# head(isoSelect(mirData, mirna="hsa-let-7a-5p", 1000))
# BiocManager::install("targetscan.Hs.eg.db")
# library("targetscan.Hs.eg.db")
# mirna_ma <- isomiRs::mirna2targetscan(c("hsa-miR-34c-5p"))
# head(mirna_ma)

#TODO miRanda,TarBase, MiRtaget2 (?)
# https://www.bioconductor.org/packages/release/bioc/html/miRNAtap.html
# BiocManager::install("miRNAtap")

#MultiQC (?) https://multiqc.info/docs/

# install.packages("sendmailR")
library(sendmailR)
from <- sprintf("sRNAflow@biomed.lu.lv","sRNAflow")
to <- sprintf(email)
subject <- paste("sRNAflow",Exp)
body <- paste("sRNAflow",Exp)
bodyWithAttachment <- list(body)
zip(file.path(ED,"fastQC.zip"),files=dir(file.path("www",Exp,"qc"),".html",full.names = T),extras="-o -j -9")
for(f in c("fastQC.zip",dir(ED,".xlsx")))
	bodyWithAttachment<-append(bodyWithAttachment,mime_part(x=file.path(ED,f),name=f))
sendmail(from,to,subject,bodyWithAttachment,control=list(smtpServer=smtpServer))



