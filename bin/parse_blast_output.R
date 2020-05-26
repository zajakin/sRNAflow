#!/usr/bin/R --no-save
source("bin/utils.R")

options(echo=TRUE)
args<-commandArgs()
if(interactive()) args[3:5] <- c(paste0("data/test"),paste0("../example-samples/Cancer_1.fa"),4)
print(args)
WD <-   args[3]
file <- args[4]
core <- args[5]

out<- paste0("blasts/",basename(file),".blast")
DV <- "nt"
blastn<- paste0("export BATCH_SIZE=50000; export BLASTDB=$HOME/data/db/blast; blastn -max_hsps 1 ")
DB <- paste("-remote -db",DV)
DB <- paste("-db",DV,"-num_threads",core)
colQuery<-  "qseqid ssciname staxid scomname sskingdom evalue bitscore qlen slen length pident mismatch qcovs stitle sseqid sstart send"
colNames<- c("read","name","taxid","nameEn","kingdom","evalue","bitscore","qlen","slen","length","pident","mismatch","qcovs","stitle","sseqid","sstart","send")
# blastopt<-paste(DB, "-evalue 1e+6 -word_size 7 -reward 2 -penalty -3 -gapopen 5 -gapextend 2")
# blastopt<-paste(DB, "-evalue 100 -perc_identity 100")
# export PATH=$HOME/bin:$HOME/conda/bin:$HOME/.local/bin:${PATH:-/usr/bin:.}; $HOME/bin/c; update_blastdb --decompress nt"
# outfmt <- "-outfmt \"6 qseqid ssciname staxid scomname evalue bitscore length pident\""
setwd(WD)
getwd()
if(!dir.exists("forKrona")) dir.create("forKrona")
if(!dir.exists("faTab")) dir.create("faTab")
if(!dir.exists("blasts")) dir.create("blasts")
if(!dir.exists("logs")) dir.create("logs")
if(!dir.exists("blasts.tmp")) dir.create("blasts.tmp")

if(length(grep(".faTab",file))==0){
    system(paste0("$HOME/conda/bin/fasta_formatter -i ",file," -t -o faTab/",basename(file),".faTab"),intern = TRUE)
    fileT<-paste0("faTab/",basename(file),".faTab")
} else fileT<-file
print(date())
if(!file.exists(paste0("forKrona/",basename(file),".forKrona.txt"))){
    bname<-paste0("blasts.tmp/",basename(file),".short.blast")
    if(!file.exists(paste0("blasts.tmp/",basename(file),".short.blast.done"))){
        system(paste0("cat ",fileT," | awk -F '\\t' '{if (length($2)<20) print \">\" $1 \"\\n\" $2}' > ",fileT,".short"),intern = TRUE)
        blastopt<-paste(DB, "-evalue 1e+6 -word_size 10 -reward 2 -penalty -3 -ungapped -perc_identity 100")
        if(file.size(paste0(fileT,".short"))>1)
            if(system(paste0(blastn,blastopt,' -outfmt \"6 ',colQuery,'\" -task blastn-short -query ',fileT,".short -out ",bname," >> logs/",basename(file),".txt 2>&1 "),intern = FALSE)==0)
                file.create(paste0("blasts.tmp/",basename(file),".short.blast.done"))
    }
    if(file.exists(bname) & file.size(bname)>1){
        short<-read.table(bname, comment.char="",quote = "", header = FALSE, sep = "\t",dec = ".", na.strings = "NA",as.is = TRUE)
        colnames(short)<-colNames
        table(short[,"qlen"]==short[,"length"])
        short<-short[short[,"qlen"]==short[,"length"],]
        write.table(short,file=bname,quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
    }
    bname<-paste0("blasts.tmp/",basename(file),".mid.blast")
    if(!file.exists(paste0("blasts.tmp/",basename(file),".mid.blast.done"))){
        system(paste0("cat ",fileT," | awk -F '\\t' '{if (length($2)>19 && length($2)<31) print \">\" $1 \"\\n\" $2}' > ",fileT,".mid"),intern = TRUE)
        blastopt<-paste(DB, "-evalue 10 -word_size 7 -reward 2 -penalty -3 -gapopen 5 -gapextend 2")
        if(file.size(paste0(fileT,".mid"))>1)
            if(system(paste0(blastn,blastopt,' -outfmt \"6 ',colQuery,'\" -task blastn-short -query ',fileT,".mid -out ",bname," >> logs/",basename(file),".txt 2>&1 "),intern = FALSE)==0)
                file.create(paste0("blasts.tmp/",basename(file),".mid.blast.done"))
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
    bname<-paste0("blasts.tmp/",basename(file),".long.blast")
    if(!file.exists(paste0("blasts.tmp/",basename(file),".long.blast.done"))){
        system(paste0("cat ",fileT," | awk -F '\\t' '{if (length($2)>30) print \">\" $1 \"\\n\" $2}' > ",fileT,".long"),intern = TRUE)
        blastopt<-paste(DB, "-evalue 0.01 -word_size 11 -reward 2 -penalty -3 -gapopen 5 -gapextend 2")
        if(file.size(paste0(fileT,".long"))>1)
            if(system(paste0(blastn,blastopt,' -outfmt \"6 ',colQuery,'\" -query ',fileT,".long -out ",bname," >> logs/",basename(file),".txt 2>&1 "),intern = FALSE)==0)
                file.create(paste0("blasts.tmp/",basename(file),".long.blast.done"))
    }
    if(file.exists(out)) file.remove(out)
    file.append(out,paste0("blasts.tmp/",basename(file),c(".short",".mid",".long"),".blast"))
}
# system(paste(blastn,blastopt,' -outfmt \"6 qseqid ssciname staxid scomname sskingdom evalue bitscore qlen slen length pident mismatch qcovs stitle\" -task blastn-short -query ',file," -out ",out," 2>&1 > ",paste0(file,".blast.log  ")),intern = TRUE)    
# bsub  -Is -M 94000 -n 4 -R "rusage[mem=94000,numcpus=4.00] span[ptile=4]" bash
    # system(paste(blastn,blastopt,' -outfmt \"6 qseqid ssciname staxid scomname sskingdom evalue bitscore qlen slen length pident mismatch qcovs stitle\" -query ',file," -out ",out),intern = TRUE )

# megablast -d nt -a $core -i $out/$i/Unmapped_$i.fasta -o $out/$i/Unmapped_$i.megablast
# blast_formatter -rid `grep -m 1 ^RID $out/$i/Unmapped_$i.blastn | awk '{print $2}'` -out $out/$i/Unmapped_$i.tab -outfmt 7 
# system("cat out.blast | sort | uniq > out_filtred.blast; $HOME/bin/parse_blast_output out_filtred.blast | sort -rn > top.txt",intern = TRUE)

# print(date())
pos<-read.table(out, comment.char="",quote = "", header = FALSE, sep = "\t",dec = ".", na.strings = "NA",as.is = TRUE)
colnames(pos)<-colNames
# head(pos)
# dim(pos)

# max(as.numeric(pos[,"evalue"]))
# pos[pos[,"evalue"]==max(pos[,"evalue"]),1:12]
# length(unique(pos[,"read"]))
# table(as.numeric(pos[!duplicated(pos[,"read"]),"qlen"])) #10-276000 11-70000 12-17700 14-1130 20-3.1   21-0.9   22-0.26
# head(pos[pos[,"qlen"]==20,])

# print(date())
species<-cbind(data.frame(unique(rbind(c("9606","Homo sapiens"),c("77133","uncultured bacterium"),pos[,3:2]))),0,0,0,0,0,NA,NA)
# species<-cbind(data.frame(cbind(unique(pos[,3]),"")),0,0,0,0)
# species[,2]<-apply(species,1,function(x) paste(pos[pos[,3]==x[1],2][1]))
rownames(species)<-species[,1]
colnames(species)<-c("id","name","pass1","pass2","count","percent","top","sstart","send")
dim(species)
# print(date())
reads<-cbind(unique(pos[,1]),0,"")
r<-reads[1,1]
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
# print(date())
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
r<-reads[1,1]
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

# reads[reads[,3]=="Esox lucius",]
# species[species[,"name"]=="Esox lucius",]
# dim(reads)
species2<-species2[order(-as.numeric(species2[,"count"])),]
head(species2)
head(reads)
species2[,"percent"]<- as.numeric(species2[,"count"])/sum(as.numeric(species2[,"count"]))*100
thr<-round(0.99*sum(as.numeric(species2[,"count"]),na.rm=TRUE))
s<-0
for(i in 1:nrow(species2)) 
    if(s>=thr){ break;
    } else s<- s+as.numeric(species2[i,"count"])
species2[1:(i-1),"top"]<- 1
species99<-species2[as.numeric(species2[,"count"])>=min(as.numeric(species2[species2[,"top"]==1,"count"])),c("id","name")]
# print(date())

fa<-read.table(fileT, comment.char="",quote = "", header = FALSE, sep = "\t",dec = ".", na.strings = "NA",as.is = TRUE)[,1]
# fa <- readLines(file)
# fa <- fa[grep("^>",fa)]
# fa <- sub("^>","",fa)
fa <- sub(" .*","",fa,perl=TRUE)
length(fa)
if(length(fa[!(fa %in% reads[,1])])>0)
    reads <- rbind(reads,cbind(fa[!(fa %in% reads[,1])],0,""))

print(date())
if(!dir.exists("species")) dir.create("species")
write.table(species2,file=paste0("species/",basename(file),".species.tsv"),quote=FALSE,row.names=TRUE,col.names=TRUE,sep="\t")
write.table(species99,file=paste0("species/",basename(file),".species99.tsv"),quote=FALSE,row.names=TRUE,col.names=TRUE,sep="\t")
if(!dir.exists("reads")) dir.create("reads")
write.table(reads,file=paste0("reads/",basename(file),".reads.tsv"),quote=FALSE,row.names=TRUE,col.names=TRUE,sep="\t")
write.table(reads[,1:2],file=paste0("forKrona/",basename(file),".forKrona.txt"),quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
if(!dir.exists("output")) dir.create("output")
system(paste0(homedir,"/conda/bin/ktImportTaxonomy forKrona/",basename(file),".forKrona.txt -o output/",basename(file),".report.htm"),intern = TRUE)
warnings()
date()
quit()

# for(i in dir("Alisa",".long.blast",recursive = TRUE)){
#     print(paste(i,length(unique(read.table(paste0("Alisa/",i), comment.char="",skip=0,quote = "", header = FALSE, sep = "\t",dec = ".", na.strings = "",as.is = TRUE)[,1])),
#         nrow(read.table(paste0("Alisa/",sub(".long.blast","",i)), comment.char="",skip=0,quote = "", header = FALSE, sep = "\t",dec = ".", na.strings = "",as.is = TRUE)) ))
# }

# Update("biomartr")
library("biomartr")
# listNCBIDatabases()
dbs<-c("refseq","genbank","ensembl","ensemblgenomes","uniprot")
# # is.genome.available(db = "refseq", species("Candidatus Carsonella ruddii PV"))
# # genomesList<-listGenomes(db = "refseq", type = "all", subset = NULL, details = FALSE)
# # length(genomesList)
# getGenome(db = "refseq", organism = "Homo sapiens")
# getGFF(db = "refseq", organism = "Homo sapiens")
# organismAttributes(organism="Homo sapiens", topic="homolog")
# getGO(organism = "Homo sapiens",genes = "GUCA2A",filters = "hgnc_symbol")


library("genomes")
genus("[Bacillus] selenitireducens")
# proks <- reports("prokaryotes.txt")
# head(proks)
species("Candidatus Carsonella ruddii PV")
species("[Bacillus] selenitireducens")
species(species99[,2])
# apply(cbind(species(species99[,2])),1,is.genome.available,db="refseq")

dir(".",".species99.tsv",recursive = TRUE)
species99<-c()
for(i in dir(".",".species99.tsv",recursive = TRUE)){
    species99<-rbind(species99,read.table(i, comment.char="",skip=0,header = TRUE, quote="",sep="\t",dec = ".", na.strings = "",as.is = TRUE))
}
species99<-species99[rownames(species99)!="77133",] # Remove "uncultered bactrium"
head(species99)
dim(species99)
species99<-unique(species99)
table(species99[,2]==species(species99[,2]))

dir.create("genomes/archive",recursive = TRUE)
dir.create("genomes/fasta",recursive = TRUE)

download.file("ftp://ftp.ensemblgenomes.org/pub/current/species.txt","genomes/species.txt")
system("sed -i '1 s/$/\tNA/' genomes/species.txt",intern = TRUE)
ensemblgenomes<-read.table("genomes/species.txt", comment.char="",skip=0,quote = "", header = TRUE, sep = "\t",dec = ".", na.strings = "",as.is = FALSE)
head(ensemblgenomes) # ,row.names = NULL 
dim(ensemblgenomes)
species(ensemblgenomes[1:10,1])

download.file("ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt","genomes/genbank.txt")
genbank<-read.table("genomes/genbank.txt", comment.char="",skip=1,quote = "", header = TRUE, sep = "\t",dec = ".", na.strings = "",as.is = TRUE)
# rownames(genbank)<-genbank[,"taxid"]
head(genbank)
dim(genbank)

download.file("ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt","genomes/refseq.txt")
refseq<-read.table("genomes/refseq.txt", comment.char="",skip=1,quote = "", header = TRUE, sep = "\t",dec = ".", na.strings = "",as.is = TRUE)
# rownames(refseq)<-refseq[,"taxid"]
head(refseq)
dim(refseq)

# refseq[1,] %in% genbank
i<-0
ext<-"_genomic.fna.gz"
file.create(paste0("genomes/genomes.fa"))
id<-"47878"
for(id in species99[,"id"]){
    print(paste(i<-i+1,id,species99[as.character(id),2],date()))
    urls<-refseq[refseq[,"taxid"]==id,"ftp_path"]
    if(length(urls)==0) urls<-refseq[refseq[,"species_taxid"]==id,"ftp_path"]
    if(length(urls)==0) urls<-genbank[genbank[,"taxid"]==id,"ftp_path"]
    if(length(urls)==0) urls<-genbank[genbank[,"species_taxid"]==id,"ftp_path"]
    if(length(urls)==0) urls<-ensemblgenomes[ensemblgenomes[,"taxonomy_id"]==id,"ftp_path"]
    dir.create(paste0("genomes/archive/",id))
    # out<-file(paste0("genomes/fasta/",id,"_",species99[as.character(id),2],".fa"),"wt")
    spname<-sub(" .*","",sub(" ","_",species99[as.character(id),2]))
    for(url in urls[1]){
        if(!file.exists(paste0("genomes/archive/",id,"/",basename(url),ext)))
            download.file(paste0(url,"/",basename(url),ext),paste0("genomes/archive/",id,"/",basename(url),ext))
        system(paste0("gunzip -c genomes/archive/",id,"/",basename(url),ext," | sed 's/^>/>",id,"_",spname,"_/g' >> genomes/genomes.fa")  )
        
    #     zz <-gzfile(paste0("genomes/archive/",id,"/",basename(url),ext),"r")
    #     # zz <-gzcon(file(paste0("genomes/archive/",id,"/",basename(urls),ext),"r"))
    #     tmp<-readLines(zz)
    #     close(zz)
    #     tmp[grep("^>",tmp)]<-sub("^>",paste0(">",id,"_",species99[as.character(id),2],"_"),tmp[grep("^>",tmp)])
    #     # cat(tmp,file=paste0(DV,".fa"),sep="\n",fill=FALSE,append=TRUE)
    #     writeLines(tmp,out,sep="\n")
    }
    # close(out)
    # file.rename(paste0(specie,ext),paste0(archiveDir,kingdom,"/",specie,ext))
    
    # ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/579/805/GCF_001579805.1_ASM157980v1
    # ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/579/805/GCF_001579805.1_ASM157980v1/GCF_001579805.1_ASM157980v1_genomic.fna.gz
}


dir.create("genomes/bowtie2")
if(file.size("genomes/genomes.fa")>(2^32-1)){ system(paste("bowtie2-build --threads ",core," --large-index genomes/genomes.fa genomes/bowtie2/genomes"))
} else system(paste("bowtie2-build --threads ",core," genomes/genomes.fa genomes/bowtie2/genomes"))

refseq[refseq[,"species_taxid"]==species99[6,"id"],"ftp_path"]
genbank[genbank[,"species_taxid"]==species99[6,"id"],"ftp_path"]

nrow(refseq[refseq[,"taxid"]==species99[6,"id"],])
nrow(refseq[refseq[,"organism_name"]==species99[6,"name"],])

table(species99[,"id"] %in% c(genbank[,"taxid"],genbank[,"species_taxid"]))

table(species99[,"id"] %in% c(refseq[,"taxid"],refseq[,"species_taxid"],genbank[,"taxid"],genbank[,"species_taxid"],ensemblgenomes[,"taxonomy_id"]))
species99[!(species99[,"id"] %in% c(refseq[,"taxid"],refseq[,"species_taxid"],genbank[,"taxid"],genbank[,"species_taxid"],ensemblgenomes[,"taxonomy_id"])),]

table(species99[,"id"] %in% refseq[,"taxid"])
table(species99[,"id"] %in% refseq[,"species_taxid"])
table(refseq[,"species_taxid"]==refseq[,"taxid"])

# _genomic.fna.gz
# _genomic.gff.gz
# _rna_from_genomic.fna.gz
# _cds_from_genomic.fna.gz
# 
# as.vector(md5sum(dir(R.home(), pattern = "^COPY", full.names = TRUE)))





