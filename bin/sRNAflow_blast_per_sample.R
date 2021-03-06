#!/usr/bin/R --no-save

print(paste(date(),"Start BLAST of subsets"))
blast_per_sample<-function(idr,re,wd,filesIn,tsize,core=4,ED){
    arg <- c("sRNAflow","--no-save",filesIn[idr,"wd"],paste0(filesIn[idr,"name"],"_random",tsize,".",re),filesIn[idr,paste0("ft",re)],core)
    options(echo=TRUE)
# if(!exists("arg")){
#     arg<-commandArgs()
#     if(interactive()) arg[3:5] <- c(paste0("www/test"),paste0("../../example-samples/Cancer_1.fa"),4)
# }
print(arg)
WD <-   arg[3]
name <- arg[4]
file <- arg[5]
core <- arg[6]

# out<- paste0("blasts/",name,".blast")
DV <- "nt"
blastn<- paste0("export BATCH_SIZE=50000; export BLASTDB=",file.path(wd,"www","db","blast"),"; blastn -max_hsps 1 -db ",DV," ")
DB <- paste("-remote")
if(file.exists(file.path(wd,"www","db","blast","db.done"))){
    DB <- paste("-num_threads",core)
    if(file.exists(file.path(wd,"www","db","meta.txids")) && blast!="nr/nt") blastn<- paste(blastn,"-taxidlist",file.path(wd,"www","db","meta.txids "))
}
colQuery<-  "qseqid ssciname staxid scomname sskingdom evalue bitscore qlen slen length pident mismatch qcovs stitle sseqid sstart send"
colNames<- c("read","name","taxid","nameEn","kingdom","evalue","bitscore","qlen","slen","length","pident","mismatch","qcovs","stitle","sseqid","sstart","send")
# blastopt<-paste(DB, "-evalue 1e+6 -word_size 7 -reward 2 -penalty -3 -gapopen 5 -gapextend 2")
# blastopt<-paste(DB, "-evalue 100 -perc_identity 100")
# export PATH=$HOME/bin:$HOME/conda/bin:$HOME/.local/bin:${PATH:-/usr/bin:.}; $HOME/bin/c; update_blastdb --decompress nt"
# outfmt <- "-outfmt \"6 qseqid ssciname staxid scomname evalue bitscore length pident\""
for(i in c("forKrona","faTab","blasts","logs")) if(!dir.exists(paste0(WD,i))) dir.create(paste0(WD,i))

setwd(WD)
getwd()
if(length(grep(".faTab",file))==0){
    fasta2tab<-" | gawk '/^>/ {printf(\"%s%s\\t\",(N>0?\"\\n\":\"\"),$1);N++;next;} {printf(\"%s\",$0);} END {printf(\"\\n\");}' > "
    system(paste0("cat ",file,fasta2tab,"faTab/",name,".faTab"),intern = TRUE)
    fileT<-paste0("faTab/",name,".faTab")
} else fileT<-file

print(date())
if(!file.exists(paste0("forKrona/",name,".forKrona.txt"))){
    bname<-paste0("blasts/",name,".short.blast")
    if(!file.exists(paste0("blasts/",name,".short.blast.done"))){
        system(paste0("cat ",fileT," | gawk -F '\\t' '{if (length($2)<20) print $1 \"\\n\" $2}' > ",fileT,".short"),intern = TRUE)
        blastopt<-paste(DB, "-evalue 1e+6 -word_size 10 -reward 2 -penalty -3 -ungapped -perc_identity 100")
        if(file.size(paste0(fileT,".short"))>1)
            if(system(paste0(blastn,blastopt,' -outfmt \"6 ',colQuery,'\" -task blastn-short -query ',fileT,".short -out ",bname," >> logs/",name,".txt 2>&1 "),intern = FALSE)==0)
                file.create(paste0("blasts/",name,".short.blast.done"))
    }
    if(file.exists(bname) & file.size(bname)>1){
        short<-read.table(bname, comment.char="",quote = "", header = FALSE, sep = "\t",dec = ".", na.strings = "NA",as.is = TRUE)
        colnames(short)<-colNames
        table(short[,"qlen"]==short[,"length"])
        short<-short[short[,"qlen"]==short[,"length"],]
        write.table(short,file=bname,quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
    }
    bname<-paste0("blasts/",name,".mid.blast")
    if(!file.exists(paste0("blasts/",name,".mid.blast.done"))){
        system(paste0("cat ",fileT," | gawk -F '\\t' '{if (length($2)>19 && length($2)<31) print $1 \"\\n\" $2}' > ",fileT,".mid"),intern = TRUE)
        blastopt<-paste(DB, "-evalue 10 -word_size 7 -reward 2 -penalty -3 -gapopen 5 -gapextend 2")
        if(file.size(paste0(fileT,".mid"))>1)
            if(system(paste0(blastn,blastopt,' -outfmt \"6 ',colQuery,'\" -task blastn-short -query ',fileT,".mid -out ",bname," >> logs/",name,".txt 2>&1 "),intern = FALSE)==0)
                file.create(paste0("blasts/",name,".mid.blast.done"))
    }
    # if(file.exists(bname) & file.size(bname)>1){
        # mid<-read.table(bname, comment.char="",quote = "", header = FALSE, sep = "\t",dec = ".", na.strings = "NA",as.is = TRUE)
        # colnames(mid)<-colNames
        # table(mid[,"qlen"]-mid[,"length"])
        # mid[mid[,"qlen"]-mid[,"length"]==11,]
        # mid[mid[,"qlen"]==19,6:12]
        # table(mid[,"mismatch"])
        # write.table(mid,file=bname,quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
    # }
    bname<-paste0("blasts/",name,".long.blast")
    if(!file.exists(paste0("blasts/",name,".long.blast.done"))){
        system(paste0("cat ",fileT," | gawk -F '\\t' '{if (length($2)>30) print $1 \"\\n\" $2}' > ",fileT,".long"),intern = TRUE)
        blastopt<-paste(DB, "-evalue 0.01 -word_size 11 -reward 2 -penalty -3 -gapopen 5 -gapextend 2")
        if(file.size(paste0(fileT,".long"))>1)
            if(system(paste0(blastn,blastopt,' -outfmt \"6 ',colQuery,'\" -query ',fileT,".long -out ",bname," >> logs/",name,".txt 2>&1 "),intern = FALSE)==0)
                file.create(paste0("blasts/",name,".long.blast.done"))
    }
    # if(file.exists(out)) file.remove(out)
    # file.append(out,paste0("blasts/",name,c(".short",".mid",".long"),".blast"))
}
# system(paste(blastn,blastopt,' -outfmt \"6 qseqid ssciname staxid scomname sskingdom evalue bitscore qlen slen length pident mismatch qcovs stitle\" -task blastn-short -query ',file," -out ",out," 2>&1 > ",paste0(file,".blast.log  ")),intern = TRUE)    
# bsub  -Is -M 94000 -n 4 -R "rusage[mem=94000,numcpus=4.00] span[ptile=4]" bash
    # system(paste(blastn,blastopt,' -outfmt \"6 qseqid ssciname staxid scomname sskingdom evalue bitscore qlen slen length pident mismatch qcovs stitle\" -query ',file," -out ",out),intern = TRUE )

# megablast -d nt -a $core -i $out/$i/Unmapped_$i.fasta -o $out/$i/Unmapped_$i.megablast
# blast_formatter -rid `grep -m 1 ^RID $out/$i/Unmapped_$i.blastn | gawk '{print $2}'` -out $out/$i/Unmapped_$i.tab -outfmt 7 
# system("cat out.blast | sort | uniq > out_filtred.blast; $HOME/bin/parse_blast_output out_filtred.blast | sort -rn > top.txt",intern = TRUE)
pos<-rbind(rep(0,length(colNames)))[-1,]
for(out in paste0("blasts/",name,c(".short",".mid",".long"),".blast"))
    pos<-rbind(pos,read.table(out, comment.char="",quote = "", header = FALSE, sep = "\t",dec = ".", na.strings = "NA",as.is = TRUE))
pos<-pos[!is.na(pos[,3]),]
colnames(pos)<-colNames

species<-cbind(data.frame(rbind(c("9606","Homo sapiens"),c("77133","uncultured bacterium"),pos[,3:2])),0,0,0,0,0,NA,NA)
species<-species[!duplicated(species[,1]),]
rownames(species)<-species[,1]
colnames(species)<-c("id","name","pass1","pass2","count","percent","top","sstart","send")
dim(species)
reads<-cbind(unique(pos[,1]),0,"")
for(r in reads[,1]){
    d<-pos[pos[,1]==r,]
    d<-d[as.numeric(d[,"evalue"])==min(as.numeric(d[,"evalue"])),]
    a<-as.character(unique(d[,"taxid"]))
    species[a,"pass1"]<- as.numeric(species[a,"pass1"])+1/length(a)
}
head(species)
species["9606","pass1"]<- as.numeric(species["9606","pass1"])+1e+6
species["77133","pass1"]<- as.numeric(species["77133","pass1"])-1e+6
species<-species[order(-as.numeric(species[,"pass1"])),]
species["9606","pass1"]<- as.numeric(species["9606","pass1"])-1e+6
species["77133","pass1"]<- as.numeric(species["77133","pass1"])+1e+6
head(species)
for(r in reads[,1]){
    b<-pos[pos[,1]==r,]
    d<-b[as.numeric(b[,"evalue"])==min(as.numeric(b[,"evalue"])) & !grepl("uncultured",b[,"name"]),]
    if(nrow(d)==0) d<-b[as.numeric(b[,"evalue"])==min(as.numeric(b[,"evalue"])),]
    a<-  as.character(species[species[,"id"] %in% d[,"taxid"],"id"])[1]
    species[a,"count"]<- as.numeric(species[a,"count"])+1
    species[a,"sstart"]<- min(as.numeric(c(d[d[,"taxid"]==a,"sstart"],species[a,"sstart"])),na.rm=TRUE)
    species[a,"send"]<- max(as.numeric(c(d[d[,"taxid"]==a,"send"],species[a,"send"])),na.rm=TRUE)
}

species2<-species[as.numeric(species[,"count"])>1 & (as.numeric(species[,"send"])-as.numeric(species[,"sstart"]))>100,]
dim(species)
dim(species2)
pos2<-pos[pos[,"taxid"] %in% species2[,"id"],]
dim(pos)
dim(pos2)
sum(species2[,"count"])
species2[,"count"]<-0
head(species2)

reads<-cbind(unique(pos2[,1]),0,"")
for(r in reads[,1]){
    d<-pos2[pos2[,1]==r,]
    d<-d[as.numeric(d[,"evalue"])==min(as.numeric(d[,"evalue"])),]
    a<-as.character(unique(d[,"taxid"]))
    species2[a,"pass2"]<- as.numeric(species2[a,"pass2"])+1/length(a)
}
head(species2)
if("9606" %in% rownames(species2)) species2["9606","pass2"]<- as.numeric(species2["9606","pass2"])+1e+6
if("77133" %in% rownames(species2)) species2["77133","pass2"]<- as.numeric(species2["77133","pass2"])-1e+6
species2<-species2[order(-as.numeric(species2[,"pass2"])),]
if("9606" %in% rownames(species2)) species2["9606","pass2"]<- as.numeric(species2["9606","pass2"])-1e+6
if("77133" %in% rownames(species2)) species2["77133","pass2"]<- as.numeric(species2["77133","pass2"])+1e+6
head(species2)

for(r in reads[,1]){
    b<-pos2[pos2[,1]==r,]
    d<-b[as.numeric(b[,"evalue"])==min(as.numeric(b[,"evalue"])) & !grepl("uncultured",b[,"name"]),]
    if(nrow(d)==0) d<-b[as.numeric(b[,"evalue"])==min(as.numeric(b[,"evalue"])),]
    a<-  as.character(species2[species2[,"id"] %in% d[,"taxid"],"id"])[1]
    species2[a,"count"]<- as.numeric(species2[a,"count"])+1
    reads[reads[,1]==r,2:3] <- as.character(species2[a,1:2])
}

# hist(species2[species2[,"count"]>1,"count"],breaks=120)
species2<-species2[order(-as.numeric(species2[,"count"])),]
species2[,"percent"]<- as.numeric(species2[,"count"])/sum(as.numeric(species2[,"count"]))*100
thr<-round(0.99*sum(as.numeric(species2[,"count"]),na.rm=TRUE))
s<-0
for(i in 1:nrow(species2)) 
    if(s>=thr){ break;
    } else s<- s+as.numeric(species2[i,"count"])
species2[1:(i-1),"top"]<- 1
species99<-species2[as.numeric(species2[,"count"])>=min(as.numeric(species2[species2[,"top"]==1,"count"])),c("id","name")]

fa<-read.table(fileT, comment.char="",quote = "", header = FALSE, sep = "\t",dec = ".", na.strings = "NA",as.is = TRUE)[,1]
fa <- sub("^>","",fa)
fa <- sub(" .*","",fa,perl=TRUE)
length(fa)
if(length(fa[!(fa %in% reads[,1])])>0)
    reads <- rbind(reads,cbind(fa[!(fa %in% reads[,1])],0,""))

if(!dir.exists("species")) dir.create("species")
write.table(species2,file=paste0("species/",name,".species.tsv"),quote=FALSE,row.names=TRUE,col.names=TRUE,sep="\t")
write.table(species99,file=paste0("species/",name,".species99.tsv"),quote=FALSE,row.names=TRUE,col.names=TRUE,sep="\t")
if(!dir.exists("reads")) dir.create("reads")
write.table(reads,file=paste0("reads/",name,".reads.tsv"),quote=FALSE,row.names=TRUE,col.names=TRUE,sep="\t")
write.table(reads[,1:2],file=paste0("forKrona/",name,".forKrona.txt"),quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
system(paste0("ktImportTaxonomy -tax ",file.path(wd,"www","db","taxonomy")," forKrona/",name,".forKrona.txt -o ",file.path(ED,"species_diagrams"),"/",name,".report.htm"),intern = TRUE)
warnings()
date()
setwd(wd)

}
