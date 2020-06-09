#!/usr/bin/R --no-save

dir(wd,".species99.tsv",recursive = TRUE)
species99<-c()
for(i in dir(wd,".species99.tsv",recursive = TRUE)){
    species99<-rbind(species99,read.table(i, comment.char="",skip=0,header = TRUE, quote="",sep="\t",dec = ".", na.strings = "",as.is = TRUE))
}
species99<-species99[rownames(species99)!="77133",] # Remove "uncultered bactrium"
print(head(species99))
print(dim(species99))
species99<-unique(species99)
#table(species99[,2]==species(species99[,2]))

