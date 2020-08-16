#!/usr/bin/R --no-save
species99<-c()
for(i in dir(ED,".species99.tsv",recursive = TRUE,full.names = TRUE)){
    species99<-rbind(species99,read.table(i, comment.char="",skip=0,header = TRUE, quote="",sep="\t",dec = ".", na.strings = "",as.is = TRUE))
}
species99<-species99[rownames(species99)!="77133",] # Remove "uncultered bactrium"
species99<-unique(species99)

