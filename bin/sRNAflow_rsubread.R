#!/usr/local/bin/Rscript --vanilla

arg<-commandArgs()
wd  <-arg[length(arg)-4]
inputFile <-arg[length(arg)-3]
gtf <-arg[length(arg)-2]
outFile  <-arg[length(arg)-1]
outSam  <-arg[length(arg)]

core   <- 78

library("Rsubread")

# genes<-featureCounts(inputFile,annot.ext=gtf,isGTFAnnotationFile=T,nthreads = core,isPairedEnd=TRUE,requireBothEndsMapped=TRUE,countMultiMappingReads=FALSE,fraction=TRUE,GTF.attrType.extra =c("gene_name","gene_biotype"))
genes<-featureCounts(inputFile,annot.ext=gtf,isGTFAnnotationFile=T,nthreads = core,isPairedEnd=FALSE,countMultiMappingReads=FALSE,fraction=TRUE,reportReads="SAM",reportReadsPath=wd,tmpDir=wd)
write.table(genes$counts[order(rownames(genes$counts)),],file = outFile, quote = FALSE, sep = "\t",col.names = FALSE)
if(outSam!="-"){
	system(paste0("grep -v 'XS:Z:Unassigned' ",inputFile,".featureCounts.sam > ",outSam," && rm ",inputFile,".featureCounts.sam"),intern = TRUE)
} else file.remove(paste0(inputFile,".featureCounts.sam"))

# d<-genes$counts
# gene_biotype<-c(genes$annotation$gene_biotype,use.names=T)
# gene_name<-c(genes$annotation$gene_name,use.names=T)
# names(gene_biotype)<-names(gene_name)<-genes$annotation$GeneID
# table(genes$annotation$GeneID==rownames(d))
# rownames(genes$stat)<-genes$stat[,1]
# # colnames(biotype)[1]<-"Status"
# table(colnames(trimm[,-1])==colnames(genes$stat[,-1]))
# stat<-rbind(trimm[,-1],.=NA,genes$stat[,-1])
# 
# biotype<-aggregate(d,list(genes$annotation$gene_biotype),sum)
# rownames(biotype)<-biotype[,1]
# proc<-cbind(round((biotype[,-1])/matrix(colSums(d),byrow=T,nrow = nrow(biotype),ncol = ncol(biotype)-1,dimnames = list(1:nrow(biotype),colnames(biotype)[-1]))*100))
# rownames(proc)<-paste0(biotype[,1],",%")
# biotype<-rbind(biotype[,-1],proc)


# htseq_opt="htseq-count -f sam -a 0 -s no --secondary-alignments score -q" #--additional-attr=gene_name
# $htseq_opt $shdir/$f.sam $DB/$DV$ext > $shdir/htseq-count_$f.txt  #  -s reverse
# Rscript --vanilla bin/rsubread.R $shdir $shdir/$f.sam $DB/$DV$ext $shdir/htseq-count_$f.txt ""

# $htseq_opt $shdir/$f.sam $shdir/tmp.gtf --samout=$shdir/priority.sam > $shdir/tmp.txt  #  -s reverse
# Rscript --vanilla bin/rsubread.R $shdir $shdir/$f.sam $shdir/tmp.gtf $shdir/tmp.txt $shdir/priority.sam

# $htseq_opt $shdir/$f.sam $DB/gtf_biotypes/$type$ext -o $shdir/tmp.sam > $shdir/htseq-nopriority_$f.$type.txt
# Rscript --vanilla bin/rsubread.R $shdir $shdir/$f.sam $DB/gtf_biotypes/$type$ext $shdir/htseq-nopriority_$f.$type.txt $shdir/tmp.sam

