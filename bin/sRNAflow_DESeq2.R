# my::Update(c("gplots","edgeR","DESeq2","biomaRt","gridExtra","ggplot2","VennDiagram"))
library(gplots)
library(corrplot)
# library(edgeR)
library(DESeq2)
library(gridExtra)
library(ggplot2)
library(VennDiagram)
tmp<-futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

options(echo=TRUE)
ED<-file.path(wd,"data",Exp)

figVen<-function(dat,lim,txt="",filexlsx){
	title<-paste0(txt," >",lim)
	vv <- list()
	for(i in 1:ncol(dat)){
		vv$tmp<-rownames(dat)[as.numeric(dat[,i])>lim]
		names(vv)[length(vv)]<-paste(colnames(dat)[i])
	}
	if(sum(lengths(vv))>0){
		aaa<-attr(venn(vv,intersections=TRUE,small=1,show.plot=FALSE),"intersections")
		delist<-as.data.frame(matrix(NA,ncol=length(aaa),nrow=max(lengths(aaa)) ))
		ccnames<-c()
		for(i in 1:length(aaa)){
			j<-length(aaa[[i]])
			if(j>0)	delist[1:j,i]<-aaa[[i]]
		}
		colnames(delist)[1:length(aaa)]<- names(aaa)
		delist[is.na(delist)]<- " "
		write.xlsx2(delist,filexlsx,sheet=txt,append=TRUE,row.names = F)
		grid.newpage()
		if(length(vv) %in% 2:3){
			if(length(vv)==2) cat.pos = c(0, 0)
			if(length(vv)==3) cat.pos = c(-40, 40, 180)
			print(grid.draw(venn.diagram(vv,fill = 2:(length(vv)+1), alpha = 0.3, filename=NULL,cex=2,cat.pos= cat.pos,
			cat.dist=0.15, main=title,cat.default.pos="outer",cat.cex=2,main.cex=2,euler.d = FALSE,scaled=FALSE)))
		} else {
			print(grid.draw(venn.diagram(vv,fill = 2:(length(vv)+1), alpha = 0.3, filename=NULL,cex=2,
			main=title,cat.default.pos="outer",cat.cex=2,main.cex=2,euler.d = FALSE,scaled=FALSE),recording = F))
		}
		wb<-openxlsx::loadWorkbook(filexlsx)
		openxlsx::insertPlot(wb,sheet=txt,width = 8, height = 6, dpi=600,startCol = 4)
		saveWorkbook(wb,filexlsx, overwrite = TRUE)
	}
}

# myEdgeR<-function(qlf,sheet,y,fileName){
# 	res<-topTags(qlf, n= nrow(y$counts) )
# 	# res<-res[as.numeric(res[[1]][,"PValue"])<1,]
# 	tmp<-(y$counts[rownames(res),]*rep(as.numeric(y$samples[,"norm.factors"]),each=nrow(res)))
# 	# colSums(tmp)
# 	# y$samples
# 	res<-cbind(rownames(res),res,NA,baseMedian=rowMedians(tmp),baseMean=rowMeans(tmp),NA,
# 		tmp,NA,y$counts[rownames(res),])
# 	write.xlsx2(as.data.frame(res),fileName,row.names=TRUE,col.names = TRUE,sheet=sheet,append=TRUE)
# }
 
makeDG<-function(httab,sel,colData,txt="test",cat1=sets[1,s],cat2=sets[2,s],filexlsx){
	dat<-httab[!(rownames(httab) %in% servRow),sel]
	dat<-dat[rowSums(dat>lim)!=0,]
	# group<-as.character(colData[sel,"nr"])
	# design <- model.matrix(~group)
	if(sum((rowSums((dat==0)+0)==0)+0)>1){
		dds <- DESeqDataSetFromMatrix(countData = dat, colData = colData[sel,], design = ~ nr)
		dds2 <- DESeq(dds)
		res <- as.data.frame(results(dds2, cooksCutoff=FALSE))
		counts<-as.matrix(counts(dds2,normalized=TRUE))[rownames(res),]
		res<-cbind(res,BH=p.adjust(res[,"pvalue"],method="BH")," "=" ",baseMedian=rowMedians(counts),"  "=" ",counts,"   "=" ",dat[rownames(res),])
		res<-res[order(res[,"pvalue"]),]
		# out<-cbind(res)
		# rownames(out)<-rownames(res)
		write.xlsx2(res,filexlsx,row.names=TRUE,col.names = TRUE,sheet=paste(txt,"DESeq2"),append=TRUE)
		out<-res[,c("log2FoldChange","pvalue","padj","baseMean")]
		colnames(out)[1:(ncol(out)-2)]<-paste(cat1,cat2,colnames(out)[1:(ncol(out)-2)])
		# res<-res[!is.na(res[,"padj"]) & res[,"padj"]<0.1,]
		# if(nrow(res)>0) write.xlsx2(as.data.frame(rbind(res,NA)),filexlsx,row.names=TRUE,col.names = TRUE,sheet=paste(txt,"DESeq2","signif"),append=TRUE)
		# res<-res[!is.na(res[,"log2FoldChange"]) & abs(res[,"log2FoldChange"])>2,]
		# if(nrow(res)>0) write.xlsx2(as.data.frame(rbind(res,NA)),filexlsx,row.names=TRUE,col.names = TRUE,sheet=paste(txt,"DESeq2","signif and Fold>2"),append=TRUE)
	} else { out<-matrix(NA,ncol=2,nrow=2); colnames(out)<-paste(paste(cat1,cat2),c("log2FoldChange","padj")) }

	# indat <- DGEList( counts = dat, group = group)
	# RLE <- calcNormFactors(indat,method="RLE")$samples[,"norm.factors"]
	# RLEgenes<-table(rowSums(dat>0)==ncol(dat))["TRUE"]
	# y <- calcNormFactors(indat,method="TMM")
	# corLibSize<-round(colSums(y$counts*rep(as.numeric(y$samples[,"norm.factors"]),each=nrow(y$counts))))
	# NFcorLibSize<-colSums(y$counts*rep(NF[sel],each=nrow(y$counts)))
	# RLEcorLibSize<-colSums(y$counts*rep(RLE,each=nrow(y$counts)))
	# foundGenes<-colSums(y$counts!=0)
	# write.xlsx2(t(cbind(y$samples,"TMM corrected libSize"=corLibSize,
	# 		"Number of genes used for RLE(DESeq2) normalisation"=RLEgenes,
	# 		"RLE NF"=RLE,"NF by mapped"=NF[sel],"NF2 by mapped to features"=NF2[sel],
	# 		"LibSize corrected by RLE"=RLEcorLibSize,"LibSize corrected by NF"=NFcorLibSize,
	# 		foundGenes=foundGenes)),filexlsx,sheet=paste(txt,"NF"),append=TRUE)
	# y <- estimateDisp(y,design)
	# myEdgeR(glmQLFTest(glmQLFit(y,design),coef=2),paste(txt,"edgeR TMM quasi-likelihood F-tests"),y,filexlsx)
	# myEdgeR(glmLRT(glmFit(y,design),coef=2),paste(txt,"edgeR TMM likelihood ratio tests"),y,filexlsx)
	# indat$samples[,"norm.factors"]<-NF[sel]
	# y<-estimateDisp(indat,design)
	# myEdgeR(glmQLFTest(glmQLFit(y,design),coef=2),paste(txt,"edgeR LibSize quasi-likelihood F-tests"),y,filexlsx)
	# indat$samples[,"norm.factors"]<-NF2[sel]
	# y<-estimateDisp(indat,design)
	# myEdgeR(glmQLFTest(glmQLFit(y,design),coef=2),paste(txt,"mappedFeat edgeR quasi-likelihood F-tests"),y,filexlsx)
	# myEdgeR(glmLRT(glmFit(y,design),coef=2),paste(txt,"edgeR LibSize likelihood ratio tests"),y,filexlsx)
	return(out)
}

stattab<-function(tt,servRow=servRow){
	stat<-rbind("Reads mapped in features"=colSums(tt[!(rownames(tt) %in% servRow),]))
	for(i in c(0,1,3,5,10,100)){
		stat<-rbind(stat,apply(tt[!(rownames(tt) %in% servRow),],2,function(x) length(x[x>i])))
		rownames(stat)[nrow(stat)]<-paste("Identifed RNA >",i,"reads")
	}
	if(length(servRow)>0) stat<-rbind(stat,"Total reads mapped"=colSums(tt))
	return(stat)
}

fileName<-file.path(ED,paste0(Exp,"_results"))
filexlsx<-paste0(fileName,".xlsx")
stat<-read.table(file.path(ED,"stats.tsv"), sep = "\t",header = TRUE)
write.xlsx2(stat[-c(21:50),],filexlsx,sheet="Quality",row.names = FALSE)
write.xlsx2(stat[21:50,],filexlsx,sheet="Catalog",row.names = FALSE,append = TRUE)

deGTF<-unique(sub(".*\\.","",sub(".txt$","",dir(ED,"_mergedFeatures.txt$",recursive = TRUE))))
prefix<-"htseq-count"
for(gr in deGTF){
	files<-file.path(wd,paste0(filesIn[,"wd"],"ShortStack/",prefix,"_",filesIn[,"name"],".",gr,".txt"))
	# DESeq2::DESeqDataSetFromHTSeqCount(files)
	httab<-cbind(read.table(files[1], sep = "\t",as.is = TRUE)[,1])
	rownames(httab)<-httab[,1]
	httab<-data.frame(httab[,-1])
	for(i in files){
		tmp<-read.table(i, header = FALSE, sep = "\t",dec = ".", na.strings = "",as.is = TRUE)
		rownames(tmp)<-tmp[,1]
		httab<-cbind(httab,as.numeric(tmp[rownames(httab),2]))
	}
	colnames(httab)<-filesIn[,"name"]
# httab<-matrix(round(runif(nrow(httab)*ncol(httab),min=-30,max = 100)),nrow=nrow(httab),ncol=ncol(httab),dimnames=list(rownames(httab),colnames(httab))); httab[httab<0]<-0
	servRow<-rownames(httab)[grep("^__",rownames(httab))]
	rr<-grep("^__",rownames(httab))

	htexp<-httab[rowSums(httab>lim)>0,]
	htexp<-rbind(htexp[!(rownames(htexp) %in% servRow),])
	print(paste(gr,grep(gr,deGTF),"/",length(deGTF),"rows =",nrow(htexp)))
	if(nrow(htexp)==0) next;
	write.xlsx2(stattab(httab,servRow),filexlsx,sheet=paste("Identifed",sub("_mergedFeatures","",gr),sep="_"),append=TRUE)

	# NF<-max(colSums(httab))/colSums(httab)
	# norm2all<-httab*max(colSums(httab))/matrix(colSums(httab),byrow=TRUE,ncol=ncol(httab),nrow=nrow(httab))
	# c2<-colSums(httab[!(rownames(httab) %in% servRow),])
	# NF2<-max(c2)/c2
	# norm2feat<-httab*max(c2)/matrix(c2,byrow=TRUE,ncol=ncol(httab),nrow=nrow(httab))
	# write.xlsx2(stattab(norm2feat,servRow),filexlsx,sheet="Normalised to features mapped reads",append=TRUE)

# unlist(unique(GroupsSel[rownames(filesIn)]))
	
	sampleTable<-cbind(Sample=colnames(httab),nr=unlist(GroupsSel[rownames(filesIn)]))
	rownames(sampleTable)<-colnames(httab)
	colData<-as.data.frame(sampleTable) # ,stringAsFactors=FALSE
	sets<-matrix(c("test","control"),nrow=2)
	for(s in 1:ncol(sets)){
		sel<-sampleTable[sampleTable[,"nr"] %in% sets[,s],"Sample"]
		out<-makeDG(httab,sel,colData,paste(sub("_mergedFeatures","",gr),sep="_"),sets[1,s],sets[2,s],filexlsx)
	}

	if(nrow(htexp)>1) figVen(htexp,lim,paste("Vennd.",sub("_mergedFeatures","",gr)),filexlsx)

	for(set in c("all", "expr.")){
		if(set=="all"){ d<-httab[!(rownames(httab) %in% servRow),] } else d<-htexp
		if(nrow(d)<3) next
		for(method in c("pearson", "spearman")){
			c<-cor(d,use="pairwise.complete.obs",method=method)
			p <- cor.mtest(d, conf.level = .95)$p
			pv<-data.frame(signif(p,2),row.names = paste0(rownames(c),"_p"),check.names = F)
			colnames(pv)<-colnames(c)
			write.xlsx2(rbind(data.frame(signif(c,2),check.names = F)," "=" ","p-values"=colnames(c),pv),
						filexlsx = filexlsx, sheet = paste(sub("_mergedFeatures","",gr),set,method), append = T)
			corrplot.mixed(c, p.mat = p, number.cex = 1, sig.level = .05,title=paste(sub("_mergedFeatures","",gr),set,method),mar=c(1,1,3,1))
			wb<-openxlsx::loadWorkbook(filexlsx)
			openxlsx::insertPlot(wb,sheet=substr(paste(sub("_mergedFeatures","",gr),set,method),0,31),width = 8, height = 6, dpi=600,startCol = 4)
			if(sum(is.na(c))==0){
				heatmap.2(c,Rowv=TRUE,Colv=TRUE, dendrogram="row", revC = T, scale="none", col=greenred(75),na.rm=TRUE, key=TRUE, density.info="none", trace="none",mar=c(8,8))
				title(paste(sub("_mergedFeatures","",gr),set,method),cex.main=0.8)
				openxlsx::insertPlot(wb,sheet=substr(paste(sub("_mergedFeatures","",gr),set,method),0,31),width = 8, height = 6, dpi=600,startCol = 15)
			}
			saveWorkbook(wb,filexlsx, overwrite = TRUE)
		}
	}
}

warnings()
date()
