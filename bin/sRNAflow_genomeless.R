#!/usr/bin/env Rscript --vanilla 
library(openxlsx)
library(EnhancedVolcano) # sudo apt install libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libproj-dev
library(PCAtools)
library(edgeR)

norm2features<-function(httab=genes[,-1],c2=colSums(httab),to=max(c2)){
	c2[c2==0]<-1
	return(httab*to/matrix(c2,byrow=TRUE,ncol=ncol(httab),nrow=nrow(httab)))
}

eger<-function(dat,colData,wb2,txt,epadj=padj,elog2FoldChange=log2FoldChange){
	keep<-rowSums((dat>lim)+0)>limS
	if(length(keep[keep])<1) return()
######################################### TODO ####################	
	# if(counts=="logical") 
	counts<-norm2features(dat[keep,])
	lib.size<-colSums(counts)
	lib.size[lib.size==0]<-1
	group<-colData[,"nr"]
	design <- model.matrix(~group)
	indat <- DGEList( counts = rbind(counts), group = group,lib.size=lib.size)
	y <- estimateDisp(indat,design)
	res <-topTags(glmQLFTest(glmQLFit(y,design),coef=2), n= nrow(y$counts) )$table #,"edgeR LibSize quasi-likelihood F-tests",y,fileName,rnk)
	# res <-topTags(glmLRT(glmFit(y,design),coef=2), n= nrow(y$counts) )$table #,"edgeR LibSize likelihood ratio tests",y,)
	if(nrow(res)==1) rownames(res)<-rownames(counts)
	colnames(res)[c(1,5)]<-c("log2FoldChange","padj")
	outDE<-cbind(res,BH=p.adjust((res[,"PValue"]),method="BH"), .=rep(" ",nrow(res)),
				 counts[rownames(res),], ..=rep(" ",nrow(res)),dat[rownames(res),])
	addWorksheet(wb = wb2, sheetName = txt, gridLines = TRUE)
	freezePane(wb2, sheet = txt, firstActiveRow = 5, firstActiveCol = 2)
	writeData(wb2, sheet = txt, paste("edgeR LibSize quasi-likelihood F-tests"),startCol=1,startRow=1)
	subtitle<-paste(txt,"significant = ",nrow(outDE[outDE[,"padj"]<epadj & abs(outDE[,"log2FoldChange"])>=elog2FoldChange,]))
	writeData(wb2, sheet = txt, subtitle,startCol=1,startRow=2)
	writeData(wb2, sheet = txt, t(colData),startCol=ncol(res)+4,colNames = F, rowNames = T,startRow=1)
	if(nrow(outDE)<500){
		fig<-tempfile()
		png(fig,width = 12, height = 12, res=150, units = "in")
		print(EnhancedVolcano(outDE[,1:6],rownames(outDE),"log2FoldChange","padj",pCutoff=epadj,
							  xlim = c(min(outDE[,"log2FoldChange"], na.rm = TRUE) - 0.5, max(outDE[,"log2FoldChange"], na.rm = TRUE) + 0.5),
							  ylim = c(0, max(-log10(outDE[,"padj"]), na.rm = TRUE) + 0.5),
							  col = c("grey30", "blue", "forestgreen", "red"),colAlpha=1,legendLabels = c("NS", expression(abs(Log[2] ~ FC) >elog2FoldChange), expression(adj. ~ "p-value < ",epadj), expression(adj. ~ "p-value < " ~ epadj ~ and ~ abs(Log[2] ~ FC) > elog2FoldChange)),
							  drawConnectors = T,maxoverlapsConnectors = 42,lengthConnectors=unit(0.01, "npc"),legendPosition = "top",pointSize = 1,subtitle=subtitle ))
		# insertPlot(wb2,sheet=txt,width = 12, height = 12, dpi=150,startCol = 20,startRow = 5)
		dev.off()
		insertImage(wb2,sheet=txt,file=fig,width = 12, height = 12, dpi=150,startCol = 20,startRow = 5)
	}
	fig<-tempfile()
	png(fig,width = 9, height = 7, res=300, units = "in")
	print(EnhancedVolcano(outDE[,1:6],NA,"log2FoldChange","padj",pCutoff=epadj,
						  xlim = c(min(outDE[,"log2FoldChange"], na.rm = TRUE) - 0.5, max(outDE[,"log2FoldChange"], na.rm = TRUE) + 0.5),
						  ylim = c(0, max(-log10(outDE[,"padj"]), na.rm = TRUE) + 0.5),
						  col = c("grey30", "blue", "forestgreen", "red"),colAlpha=1,legendLabels = c("NS", expression(abs(Log[2] ~ FC) ~ "> 1"), expression("adj. p-value < 0.05"), expression(abs(Log[2] ~ FC) ~ "> 1 and adj. p-value < 0.05")),
						  drawConnectors = F,lengthConnectors=unit(0.01, "npc"),legendPosition = "top",pointSize = 1,subtitle=subtitle))
	# insertPlot(wb2,sheet=txt,width = 9, height = 7, dpi=300,startCol = 10,startRow = 5)
	dev.off()
	insertImage(wb2,sheet=txt,file=fig, width = 9, height = 7, dpi=300,startCol = 10,startRow = 5)
	if(nrow(indat)>1){
		# plotMDS(indat,main = "Distances on the plot approximate the typical log2 fold changes between the samples\ntop genes separately for each pairwise comparison between the samples")
		# insertPlot(wb2,sheet=txt,width = 9, height = 7, dpi=150,startCol = 20,startRow = 63)
		# plotMDS(indat, gene.selection="common",main = "Distances on the plot approximate the typical log2 fold changes between the samples\nPCA like: the same genes for all comparisons")
		# insertPlot(wb2,sheet=txt,width = 9, height = 7, dpi=150,startCol = 20,startRow = 98)
		p <- pca(cpm(indat), metadata = colData) ## , removeVar = 0.1 -- removing the lower 10% of variables based on variance
		fig<-tempfile()
		png(fig, width = 9, height = 7, res=150, units = "in")
		print(screeplot(p, axisLabSize = 18, titleLabSize = 22))
		# insertPlot(wb2,sheet=txt, width = 9, height = 7, dpi=150,startCol = 10,startRow = 75)
		dev.off()
		insertImage(wb2,sheet=txt,file=fig, width = 9, height = 7, dpi=150,startCol = 10,startRow = 75)
		fig<-tempfile()
		png(fig, width = 9, height = 7, res=300, units = "in")
		print(biplot(p, colby = 'nr', colLegendTitle = txt, encircle = TRUE, encircleFill = TRUE, hline = 0, vline = c(-25, 0, 25),
					 legendPosition = 'top', legendLabSize = 16, legendIconSize = 8.0, showLoadings = TRUE, sizeLoadingsNames = 5))
		# insertPlot(wb2,sheet=txt,width = 9, height = 7, dpi=300,startCol = 10,startRow = 40)
		dev.off()
		insertImage(wb2,sheet=txt,file=fig, width = 9, height = 7, dpi=300,startCol = 10,startRow = 40)
	}
	outDE<-outDE[outDE[,"padj"]<epadj & abs(outDE[,"log2FoldChange"])>=elog2FoldChange,]
	writeData(wb2, sheet = txt, outDE[1:min(nrow(outDE),3000),], colNames = T, rowNames = T,startRow=4) # [order(outDE[,colnames(outDE)[grep ("t.test",colnames(outDE))[1]]]),]
	return(outDE)
}

myphylogenicTree<-function(dl,tc,methodCl="ClustalW",rowLimit=500,txt="",wb2=wb,sheet=sheet){ # "ClustalW", "Muscle"
	library(msa)
	library(seqinr)
	library(ape)
	dl<-dl[!is.na(dl[,1]),]
	fn<-file.path(tc,paste0(txt,"_",min(rowLimit,nrow(dl))))
	if(!file.exists(paste0(fn,".fa"))){
		tmp<-data.frame(Id=NA,dl[1:min(rowLimit,nrow(dl)),],check.names = F)
		DB <- paste("-remote")
		if(file.exists(file.path(wd,"www","db","blast","db.done"))){
			DB <- paste("-num_threads",core)
			if(file.exists(file.path(wd,"www","db","meta.txids")) && blast!="nr/nt") DB<- paste(DB,"-taxidlist",file.path(wd,"www","db","meta.txids "))
		}
		# cat(c(rbind(paste0(">",rownames(tmp)),rownames(tmp))),file=paste0(fn,"_4blast.fa"), sep="\n")
		colQuery<-  "qseqid ssciname staxid scomname sskingdom evalue bitscore qlen slen length pident mismatch qcovs stitle sseqid sstart send"
		# colNames<- c("read","name","taxid","nameEn","kingdom","evalue","bitscore","qlen","slen","length","pident","mismatch","qcovs","stitle","sseqid","sstart","send")
		blastn<-paste0("export BATCH_SIZE=50000; export BLASTDB=",file.path(wd,"www","db","blast"),"; blastn -max_hsps 1 -db nt ",DB)
		
		cat(c(rbind(paste0(">",rownames(tmp)[nchar(rownames(tmp))<20],"_"),rownames(tmp)[nchar(rownames(tmp))<20])),file=paste0(fn,"_short.fa"), sep="\n")
		blastopt<-" -evalue 1e+6 -word_size 10 -reward 2 -penalty -3 -ungapped -perc_identity 100"
		if(file.size(paste0(fn,"_short.fa"))>1)
			system(paste0(blastn,blastopt,' -outfmt \"6 ',colQuery,'\" -query ',paste0(fn,"_short.fa")," -out ",paste0(fn,"_short.tsv")," > ",paste0(fn,"_blast.log")," 2>&1 "),intern = FALSE)
		
		cat(c(rbind(paste0(">",rownames(tmp)[nchar(rownames(tmp))>19 & nchar(rownames(tmp))<31],"_"),rownames(tmp)[nchar(rownames(tmp))>19 & nchar(rownames(tmp))<31])),file=paste0(fn,"_middle.fa"), sep="\n")
		blastopt<-" -evalue 10 -word_size 7 -reward 2 -penalty -3 -gapopen 5 -gapextend 2"
		if(file.size(paste0(fn,"_middle.fa"))>1)
			system(paste0(blastn,blastopt,' -outfmt \"6 ',colQuery,'\" -query ',paste0(fn,"_middle.fa")," -out ",paste0(fn,"_middle.tsv")," >> ",paste0(fn,"_blast.log")," 2>&1 "),intern = FALSE)
		
		cat(c(rbind(paste0(">",rownames(tmp)[nchar(rownames(tmp))>30]),rownames(tmp)[nchar(rownames(tmp))>30],"_")),file=paste0(fn,"_long.fa"), sep="\n")
		blastopt<-" -evalue 0.01 -word_size 11 -reward 2 -penalty -3 -gapopen 5 -gapextend 2"
		if(file.size(paste0(fn,"_long.fa"))>1)
			system(paste0(blastn,blastopt,' -outfmt \"6 ',colQuery,'\" -query ',paste0(fn,"_long.fa")," -out ",paste0(fn,"_long.tsv")," >> ",paste0(fn,"_blast.log")," 2>&1 "),intern = FALSE)
		
		# system(paste0(opt,"echo ",i," | blastn -max_hsps 1 -db nt -num_threads ",core,' -outfmt \"6 ',colQuery,'\" -query ',paste0(fn,"_4blast.fa")," -out ",paste0(fn,"_blast.tsv")," >> ",paste0(fn,"_blast.log")," 2>&1 "),intern = FALSE)
		for(i in c("short","middle","long")) if(!file.exists(paste0(fn,"_",i,".tsv")) | file.size(paste0(fn,"_",i,".tsv"))==0) cat(paste(rep("-",length(strsplit(colQuery," ")[[1]])),collapse = "\t"),file=paste0(fn,"_",i,".tsv"), sep="\n")
		bl<-rep("-",length(strsplit(colQuery," ")[[1]]))
		for(i in c("short","middle","long")) bl<-rbind(bl,read.table(paste0(fn,"_",i,".tsv"),sep="\t",quote = "",comment.char = "") )
		bl[,1]<-sub("_$","",bl[,1])
		head(bl)
		tmp[,"Id"]<-unlist(lapply(rownames(tmp),function(x){ aaa<-bl[bl[,1]==x[1],]; paste0(aaa[1,2]," ",aaa[1,16]," [",aaa[1,5],"] ",aaa[1,11],"%") }))
		tmp[tmp[,"Id"]=="NA NA [NA] NA%","Id"]<-"No hits"
		cat(c(rbind(paste0(">",rownames(tmp)," ",round(tmp[,"log2FoldChange"],1)," ",signif(tmp[,"padj"],2)," ",apply(tmp[,(grep("..",colnames(tmp),fixed=T)+1):ncol(tmp)],1,sum)," ",tmp[,"Id"]),rownames(tmp))),file=paste0(fn,".fa"), sep="\n")
	}
	if(!file.exists(paste0(fn,"_",methodCl,".RData"))){
		aaa<-msaConvert(msa(readDNAStringSet(paste0(fn,".fa")),method=methodCl), type="seqinr::alignment")
		dist<-dist.alignment(aaa, "identity")
		njs<-njs(dist)
		save(aaa,njs,dist,file=paste0(fn,"_",methodCl,".RData"))
	} else  load(paste0(fn,"_",methodCl,".RData"))
	# par(family="mono",oma=c(0,0,0,0),mar=c(1,1,1,1))
	colourtips<-rep("black",length(njs$tip.label))
	colourtips[grep("Bacteria",njs$tip.label,ignore.case = T)]<-"darkgreen"
	colourtips[grep("Eukaryota",njs$tip.label,ignore.case = T)]<-"red"
	colourtips[grep("fung",njs$tip.label,ignore.case = T)]<-"darkmagenta"
	colourtips[grep("[ATGCN] -[0-9]",njs$tip.label)]<-"blue"

	fig<-tempfile()
	png(fig, width = 57, height = 4+length(njs$tip.label)*0.21, res=75, units = "in",family="mono")
	plot.phylo(njs, main=methodCl,no.margin = T,align.tip.label=T,tip.color=colourtips,use.edge.length=T)
	# insertPlot(wb2,sheet=sheet,width = 57, height = 4+length(njs$tip.label)*0.21, dpi=75,startCol = 6,startRow = 2)
	dev.off()
	insertImage(wb2,sheet=sheet,file=fig, width = 57, height = 4+length(njs$tip.label)*0.21, dpi=75,startCol = 6,startRow = 2)
	# par(family="",oma=c(0,0,0,0),mar=c(5.1,4.1,4.1,2.1))
	return(setNames(aaa$seq,aaa$nam)[njs$tip.label])
}

genomless<-function(ED,colData){
	wb<-createWorkbook()
	tc<-file.path(ED,paste0("genomeless",tax))
	if(!dir.exists(tc)) dir.create(tc,recursive = T)
	if(!file.exists(file.path(tc,"gless.RData"))){
		ext<-paste0(".reads.gz")
		if(strategy=="successively"){
			foreach(f=colData$Sample) %dopar% system(paste0("sed -n '2~4p' ",file.path(ED,f,""),"Unmapped_2main_",f,".fq | sort | uniq -c | pigz -9cf > ",file.path(ED,f,f),tax,".reads.gz"),intern = T)
			ext<-paste0(tax,".reads.gz")
		}
		keep<-c()
		for(f in colData$Sample){
			sk<-data.frame(read.table(paste0(file.path(ED,f,f),ext)))
			if(length(sk)==2) keep<-c(keep,sk[sk[,1]>lim,2])
		}
		tmp<-table(keep)
		keep<-names(tmp[tmp>limS])
		if(exists("tabgl")) rm("tabgl")
		tabgl<-matrix(0,nrow = length(keep), ncol = length(colData$Sample) , dimnames = list(keep,colData$Sample))
		for(f in colData$Sample){
			sk<-data.frame(read.table(paste0(file.path(ED,f,f),ext)),row.names = 2)
			overlap<-keep[keep %in% rownames(sk)]
			if(length(sk)==1) tabgl[overlap,f]<-sk[overlap,1]
		}
		tabgl[is.na(tabgl)]<-0
		# tabgl<-tabgl[rowSums(tabgl>lim)>limS,]
		save(tabgl,file=file.path(tc,"gless.RData"))
	} else load(file.path(tc,"gless.RData"))
	# paste("How many hits present in all samples (NA - means not present) : ",table(rowSums(tabgl==0)>0)["FALSE"])
	nfGL<-10^6/colSums(tabgl)
	gless<-eger(dat=tabgl,colData=colData,wb2=wb,txt=paste(Exp,"genomeless"))
	# save(gless,file=file.path(tc,"gless2.RData")) # load(file.path(tc,"gless2.RData"))
	if(nrow(gless)>0){
		sheet<-methodCl<-"ClustalW"
		addWorksheet(wb,sheet)
		modifyBaseFont(wb, fontSize = 11, fontColour = "black", fontName = "Courier New")
		setColWidths(wb,sheet,cols=1:4,widths = "auto")
		rowLimit<-2000
		clusters<-data.frame(myphylogenicTree(dl=gless[c(rownames(gless)[gless[,"log2FoldChange"]>0][1:min(rowLimit/2,nrow(gless[gless[,"log2FoldChange"]>0,]))],rownames(gless)[gless[,"log2FoldChange"]<0][1:min(rowLimit/2,nrow(gless[gless[,"log2FoldChange"]<0,]))]),],tc,methodCl,rowLimit=rowLimit,txt = paste0("Clusters"),wb,sheet=sheet))
		clusters<-data.frame(ClustalW=clusters[,1], Seq=sub(" .*","",rownames(clusters)), BLAST=sub(".*%%%","",sub(" ","%%%",rownames(clusters))))
		writeDataTable(wb,sheet,clusters,rowNames = T, colNames = T)
	}
	# tc<-file.path(ED,paste0("genomeless",tax),paste0(cancerTypes[project,1],"_",comp))
	# if(!dir.exists(tc)) dir.create(tc,recursive = T)
	DV<-file.path(tc,"genome",paste0(Exp,"_Genomeless.fa"))
	if(!dir.exists(file.path(tc,"genome"))) dir.create(file.path(tc,"genome"))
	if(file.exists(DV)){
		if(!dir.exists(file.path(tc,"mapped"))) dir.create(file.path(tc,"mapped"))
		if(!file.exists(sub(".fa$",".1.bt2",DV))) system(paste0("bowtie2-build --threads ",core," ",DV," ",sub(".fa$","",DV)," > ",sub(".fa$",".log",DV)))
		bowtie2opt<-paste("--time --end-to-end -p ",core," --mm --no-unal --very-sensitive --no-1mm-upfront --score-min L,-1.15,-0.24")

		foreach(f=colData$Sample) %dopar% if(!file.exists(paste0(file.path(tc,"mapped",f),".txt")))
			system(paste0("bowtie2 ",bowtie2opt," -x ",sub(".fa$","",DV)," -U ",ED,"/",f,"/Unmapped_2main_",f,".fq -S ",file.path(tc,"mapped",f),".sam > ",file.path(tc,"mapped",f),".log 2>&1 && ",
						  "cat ",file.path(tc,"mapped",f),".sam | grep -v -F -f ",ED,"/genomeless",tax,"/exclude | samtools view -uhS -F4 | samtools sort -@ ",core," - -o ",file.path(tc,"mapped",f),".bam >> ",file.path(tc,"mapped",f),".log 2>&1 && ", # -w
						  # "~/igv/IGV/igvtools index ",file.path(tc,f),".sam >> ",file.path(tc,f),".log 2>&1 && ",
						  "samtools index ",file.path(tc,"mapped",f),".bam  >> ",file.path(tc,"mapped",f),".log 2>&1 && ",
						  "samtools idxstats ",file.path(tc,"mapped",f),".bam  > ",file.path(tc,"mapped",f),".txt"))
		rows<-sub("^>","",system(paste0("grep '>' ",DV),intern = T))
		tabGG<-data.frame(row.names = rows)
		for(f in colData$Sample){
			tmp<-read.table(paste0(file.path(tc,"mapped",f),".txt"),row.names = 1)
			tmp<-tmp[rownames(tmp)!="*",]
			tabGG[[f]]<-tmp[rows,2]
		}
		tabGG[is.na(tabGG)]<-0
		gg<-eger(dat=tabGG[rowSums(tabGG>0)>1,],colData=colData,wb2=wb,txt=paste0("GL_selected"),
				 counts=tabGG[rowSums(tabGG>0)>1,]*matrix(nfGL[colnames(tabGG)],byrow=TRUE,ncol=ncol(tabGG),nrow=nrow(tabGG[rowSums(tabGG>0)>1,])))
		# toIGV_GL(gg,comp,ED,colnames(tabgl))
		if(!file.exists(sub(".fa$",".html",DV))){
			DB <- paste("-remote")
			if(file.exists(file.path(wd,"www","db","blast","db.done"))){
				DB <- paste("-num_threads",core)
				if(file.exists(file.path(wd,"www","db","meta.txids")) && blast!="nr/nt") DB<- paste(DB,"-taxidlist",file.path(wd,"www","db","meta.txids "))
			}
			blastn<- paste0("export BATCH_SIZE=50000; export BLASTDB=",file.path(wd,"www","db","blast"),"; blastn -max_hsps 1 -db nt ",DB)
			system(paste0(blastn," -query ",DV," -html -out ",sub(".fa$",".html",DV)),intern = FALSE)
		}
		rowSums(tabGG)[order(rowSums(tabGG))]
		tabGG<-tabGG[order(rowSums(tabGG)),]
		
		tmp<-gg[,(grep("^\\.$",colnames(gg))+1):(grep("^\\.\\.$",colnames(gg))-1)]
		tmp<-apply(tmp,1,function(x) round(as.numeric(performance(prediction(as.numeric(x),colData[colnames(tmp),"nr"]), measure = "auc", x.measure = "cutoff")@y.values),digits=3))
		
		sheet<-paste0("AUC")
		addWorksheet(wb,sheet)
		writeDataTable(wb,sheet,data.frame(Sequence=sub("_.*","",rownames(gg)),Specie=sub(paste0(".*_vs_[A-Za-z]+_"),"",rownames(gg)),gg[,"log2FoldChange"],gg[,"padj"],AUC=tmp,Length=nchar(sub("_.*","",rownames(gg))) ),rowNames = FALSE, colNames = TRUE)
	}
	saveWorkbook(wb,file.path(ED,paste0(Exp,"_genomeless.xlsx")),overwrite = TRUE)
}

colData<-data.frame(Sample=filesIn[,"name"],nr=filesIn[,"gr"],row.names = filesIn[,"name"])
colData<-colData[colData$nr %in% c("test","control"),]
# colData$nr[colData$nr=="test"]   <-1
# colData$nr[colData$nr=="control"]<-0
genomless(ED,colData)

