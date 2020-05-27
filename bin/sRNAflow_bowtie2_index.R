#!/usr/bin/env Rscript --vanilla 
options(echo=TRUE)
args<-commandArgs()
if(interactive()) args[3:4] <- c("data/test/genomes/genomes.fa",4)
print(args)
fa <-   args[3]
core<-4 
if(args[4] != "") core<-args[4]
if(file.size(fa)>(2^32-1)){ system(paste("bowtie2-build --threads",core," --large-index ",fa,sub(".fa","",fa)))
} else system(paste("bowtie2-build --threads",core,fa,sub(".fa","",fa)))
