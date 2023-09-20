library(gplots)
library(corrplot)
# library(edgeR)
library(DESeq2)
library(gridExtra)
library(ggplot2)
library(VennDiagram)
tmp<-futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

options(echo=TRUE)
setwd(file.path(ED,"species_diagrams"))
print(paste(date(),"Compress species diagrams"))
zip(file.path(ED,paste0(Exp,"_sRNAflow_diagrams.zip")),files=dir(".",include.dirs = TRUE),flags="-roq9")
con<-file(file.path("..",paste0(Exp,"_species_diagrams.html")),"wt")
cat(paste0('<html><head><link rel="stylesheet" type="text/css" href="/shared/shiny.css"/></head><body><table width="100%" border="2">'),file = con)
index<-dir(".",".htm$")
index<-index[!grepl(paste0("_random",tsize,"."),index)]
for(i in index){
	cat(paste0("<tr align=center><td><a href=species_diagrams/",i,">",i,"</a></td>"),file = con)
	for(j in 1:Rep) cat(paste0("<td><a href=species_diagrams/",sub(".report.htm",paste0("_random",tsize,".",j,".report.htm"),i),">",
							   sub(".report.htm",paste0("_random",tsize,".",j,".report.htm"),i),"</a></td>"),file = con)
	cat(paste("</tr>"),file = con)
}
cat(paste("</table></body></html>"),file = con)
close(con)
setwd(wd)

print(paste(date(),"isomiR-SEA"))
tabs<-c(c("tag_unique","ambigue","unique"),paste(c("tag","ambigue","unique"),"ambigue_selected",sep = "_")) #"tag_ambigue",
if(!dir.exists(file.path(ED,"isomiR-SEA"))) dir.create(file.path(ED,"isomiR-SEA"),recursive = TRUE, mode = "0777")
for(sample in filesIn[,"name"]){
	p<- file.path(ED,sample,"isomiR-SEA")
	filexlsx<- file.path(ED,"isomiR-SEA",paste0(sample,"_isomiR-SEA.xlsx"))
	tab<-read.table(file.path(p,"summary.txt"), sep = " ",header = F)[,-2]
	write.xlsx2(tab,filexlsx,sheet="Summary",row.names = FALSE,col.names = FALSE)
	for(t in tabs){
		infile<-file.path(p,paste0("out_result_mature_",t,".txt"))
		if(length(readLines(infile,n=2))==2){
			tab<-read.table(infile, sep = "\t",header = F,skip = 1,comment.char = "",fill = T)
			colnames(tab)<-read.table(infile,header = FALSE,nrows = 1,skip = 0,comment.char = "")
			write.xlsx2(tab,filexlsx,sheet=t,row.names = FALSE,append = TRUE)
		}
	}	
}
zip(file.path(ED,paste0(Exp,"_isomiR-SEA.zip")),files=dir(file.path(ED,"isomiR-SEA"),"_isomiR-SEA.xlsx",full.names = T,recursive = TRUE),flags="-joq9")

figVen<-function(dat,lim,txt="",filexlsx){
	title<-paste0(txt," >",lim)
	vv <- list()
	if(ncol(dat)>5 || ncol(dat)<2) return()
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
		openxlsx::insertPlot(wb,sheet=txt,width = 8, height = 6, dpi=150,startCol = 4)
		saveWorkbook(wb,filexlsx, overwrite = TRUE)
	}
}

tomieaa<-function(mirnas,ssp="hsa",analysis="ora"){
	mirnas<-mirnas[grep(paste0("^",ssp,"-"),mirnas)]
	mature<-mirnas[grep("p$",mirnas)]
	precursors<-mirnas[!grepl("p$",mirnas)]
	if(length(mature)>0) precursors<-append(precursors,strsplit(system(paste("mieaa to_precursor -m ",mature," --unique"),intern = TRUE),"'")[[1]][c(FALSE,TRUE )])
	system(paste("mieaa",analysis,ssp,"--precursors -m ",paste(precursors,collapse=",")," -c HMDD MNDR -o ",file.path(ED,"mieaa.csv")),intern = TRUE)
	df<-"Not found"
	if(file.size(file.path(ED,"mieaa.csv"))>0) df<-read.csv(file.path(ED,"mieaa.csv"))
	file.remove(file.path(ED,"mieaa.csv"))
	return(df)
}

myGO<-function(myids, minimum.population=5,deUp=c(),deDown=c()){
	library(org.Hs.eg.db)
	library(GOstats)
	keyunis <- keys(org.Hs.eg.db, keytype="GOALL")
	output <- data.frame(total=integer(0), expected=numeric(0), observed=integer(0), p.value=numeric(0),adj_pval=numeric(0), 
						 OddsRatio=character(0),Term=character(0), ontology=character(0), All_genes=character(0), count_DE_genes=integer(0),
						 DE_genes=character(0),count_Up=integer(0), count_Down=integer(0),genes_Up=character(0), genes_Down=character(0))
	oName<-matrix(c('BP','biological process', 'MF','molecular function','CC','cellular component'),nrow=2)
	for(i in 1:3){
		params <- new('GOHyperGParams', geneIds=myids, annotation="org.Hs.eg.db", ontology=oName[1,i],
					  pvalueCutoff=0.05, conditional=TRUE, testDirection="over")
		go <- hyperGTest(params)
		# gn2go <- select(org.Hs.eg.db, as.character(geneIds(go)), "GOALL")
		gn <- select(org.Hs.eg.db, unique(as.character(geneIds(go))), "SYMBOL")
		rownames(gn)<-gn[,1]
		go.table <- summary(go, pvalue=2)
		# ggo <- select(org.Hs.eg.db, rownames(go.table), "GOALL")
		ggo <- select(org.Hs.eg.db, rownames(go.table), keys=keyunis, columns = c("GOALL","SYMBOL","ENTREZID"), keytype="GOALL")
		# columns(org.Hs.eg.db)
		if (nrow(go.table)>0) {
			go.table$Pvalue <- p.adjust(go.table$Pvalue, method="none")
			go.table$adj_pval <- p.adjust(go.table$Pvalue, method="fdr")
			go.table <- go.table[go.table$Pvalue <= 0.05 & go.table$Size >= minimum.population,]
			if (nrow(go.table)>0) {
				rownames(go.table) <- go.table[,1]
				go.table <- go.table[,c(6,4,5,2,8,3,7)]
				go.table$ontology <- oName[1,i]
				for(g in 1:nrow(go.table)){
					sid<-which(ggo[,"GOALL"]==rownames(go.table)[g])
					names<-unique(ggo[sid,"SYMBOL"])
					# names<-gn[gn2go[sid,1],2]
					go.table[g,"All_genes"]<-paste(names[order(names)],collapse=",")
					names<-gn[as.character(myids[myids %in% ggo[sid,"ENTREZID"]]),2]
					up<-gn[as.character(deUp[deUp %in% ggo[sid,"ENTREZID"]]),2]
					down<-gn[as.character(deDown[deDown %in% ggo[sid,"ENTREZID"]]),2]
					# names<-gn[as.character(myids[myids %in% gn[gn2go[sid,1],1]]),2]
					go.table[g,"count_DE_genes"]<-length(names)
					go.table[g,"DE_genes"]<-paste(names,collapse=",")
					go.table[g,"count_Up"]<-  length(up)
					go.table[g,"count_Down"]<-length(down)
					go.table[g,"genes_Up"]<-  paste(up,collapse=",")
					go.table[g,"genes_Down"]<-paste(down,collapse=",")
				}
				colnames(go.table) <- colnames(output)
				output <- rbind(output, go.table)
			}
		}
	}
	return(output[order(output$p.value),])
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
 
makeDG<-function(httab,sel,colData,txt="test",cat1=sets[1,s],cat2=sets[2,s],filexlsx,servRow=c()){
	dat<-httab[!(rownames(httab) %in% servRow),sel]
	dat<-as.matrix(dat[rowSums(dat>lim)!=0,])
	mode(dat)<-"integer"
	# group<-as.character(colData[sel,"nr"])
	# design <- model.matrix(~group)
	colData[,2]<-as.factor(colData[,2])
	if(sum((rowSums((dat==0)+0)==0)+0)>1){
		dds <- DESeqDataSetFromMatrix(countData = dat, colData = colData[sel,], design = ~ nr)
		err<-try(dds2 <- DESeq(dds,quiet = T))
		if(typeof(err)!="S4"){
			write.xlsx2(gsub("\n"," ",as.character(err)),filexlsx,row.names=TRUE,col.names = TRUE,sheet=paste(txt,"DESeq2"),append=TRUE)
			return(rbind(list(baseMean=1,log2FoldChange=1,lfcSE=1,stat=1,pvalue=1,padj=1))[-1,])
		}
		res <- as.data.frame(results(dds2,contrast=c("nr","test","control"),lfcThreshold=log2FoldChange, alpha=padj, cooksCutoff=TRUE))
		res<-res[order(res[,"pvalue"]),]
		out<-res
		res<-res[!is.na(res[,"padj"]) & res[,"padj"]<padj,]
		if(nrow(res)>0){
			counts<-rbind(as.matrix(counts(dds2,normalized=TRUE))[rownames(res),])
			res<-cbind(rbind(res)," "=" ",baseMedian=rowMedians(counts),"  "=" ",counts,"   "=" ",Raw=rbind(dat[rownames(res),]))
			write.xlsx2(as.data.frame(rbind(res,NA)),filexlsx,row.names=TRUE,col.names = TRUE,sheet=paste(txt,"DESeq2"),append=TRUE)
			#TODO Make working for mouse etc...  ####
			if(txt == "all_miRNA") try(write.xlsx2(tomieaa(sub("\\[[+-.]\\]","",unlist(strsplit(sub(".*_mergedFeatures_","",rownames(res)),"/"))),analysis="ora"),
												   filexlsx,row.names=TRUE,col.names = TRUE,sheet=paste(txt,"MIEAA_ORA"),append=TRUE))
			if(txt == "all_miRNA") try(write.xlsx2(tomieaa(sub("\\[[+-.]\\]","",unlist(strsplit(sub(".*_mergedFeatures_","",rownames(res)),"/"))),analysis="gsea"),
												   filexlsx,row.names=TRUE,col.names = TRUE,sheet=paste(txt,"MIEAA_GSEA"),append=TRUE))
			#TODO Diff expression DESeq2  GOstats ####
			# if(txt %in% c("Ensembl_genes","protein_coding","all_lncRNA","all_miRNA")){
			# # 	#  https://www.bioconductor.org/packages/release/bioc/manuals/GOstats/man/GOstats.pdf
			# 	myGO(myids=sub("\\[[+-.]\\]","",unlist(strsplit(sub(".*_mergedFeatures_","",rownames(res)),"/"))), minimum.population=5,deUp=c(),deDown=c())
			# }
		} else write.xlsx2(c("No matching results found"),filexlsx,row.names=TRUE,col.names = TRUE,sheet=paste(txt,"DESeq2"),append=TRUE)
	} else write.xlsx2(c("Normalisation not possible: No found features presented in all samples"),filexlsx,row.names=TRUE,col.names = TRUE,sheet=paste(txt,"DESeq2"),append=TRUE)

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

print(paste(date(),"Differencial expressions tables"))
fileName<-file.path(ED,paste0(Exp,"_results"))
filexlsx<-paste0(fileName,".xlsx")
value<-c(Exp=Exp,specie=specie,tsize=tsize,Rep=Rep,blast=blast,ad3=ad3,ad5=ad5,sizerange =sizerange,lim =lim,log2FoldChange=log2FoldChange,padj =padj,email =email,smtpServer=smtpServer)
value<-rbind(cbind(variable=names(value),value=value," "=" ","  "=" ")," ",c("file","size","date","group"),cbind(FilesIn,group=GroupsSel[FilesIn[,"file"]]))
write.xlsx2(value,filexlsx,sheet="Settings",row.names = FALSE)
stat<-read.table(file.path(ED,paste0("stats",tax,".tsv")), sep = "\t",header = TRUE)
enrow<-grep("all_Ensembl",stat[,1])
write.xlsx2(stat[1:enrow,],filexlsx,sheet="Quality",row.names = FALSE,append = TRUE)
write.xlsx2(stat[(enrow+1):nrow(stat),],filexlsx,sheet="Catalog",row.names = FALSE,append = TRUE)

endrow<- grep("Ensembl_genes_mergedFeatures",stat[,1])-1
samples <- rep(colnames(stat)[-1],each=(endrow-enrow))
RNA_types <- rep(stat[(enrow+1):endrow,1], ncol(stat)-1)
frequency <- as.numeric(unlist(stat[(enrow+1):endrow,-1]))
data <- data.frame(samples,RNA_types,frequency)
png(paste0(fileName,".png"),width = 6+ncol(stat)/4, height = 8, res=300, units = "in")
print(ggplot(data, aes(fill=RNA_types, y=frequency, x=samples)) + geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()
wb<-openxlsx::loadWorkbook(filexlsx)
# openxlsx::insertPlot(wb,sheet="Catalog",width = 6+ncol(stat)/4, height = 8, dpi=300,startCol = 8)
openxlsx::insertImage(wb,sheet="Catalog",file=paste0(fileName,".png"),width = 6+ncol(stat)/4, height = 8, dpi=300,startCol = 8)
saveWorkbook(wb,filexlsx, overwrite = TRUE)
file.remove(paste0(fileName,".png"))
# deGTF<-unique(sub(".*\\.","",sub(".txt$","",dir(ED,"_mergedFeatures.txt$",recursive = TRUE))))
deGTF<-c("all_miRNA_mergedFeatures","GtRNAdb_mergedFeatures","rRNA_mergedFeatures","protein_coding_mergedFeatures","processed_pseudogene_mergedFeatures","snRNA_mergedFeatures","snoRNA_mergedFeatures","all_MT_mergedFeatures","all_piRNA_mergedFeatures","all_lncRNA_mergedFeatures","vault_RNA_mergedFeatures","misc_RNA_mergedFeatures","Other_types_mergedFeatures","RepeatMasker_tRNA_mergedFeatures","RepeatMasker_rRNA_mergedFeatures")

prefix<-"htseq-count"
for(gr in deGTF){
	files<-paste0(filesIn[,"wd"],"ShortStack/",prefix,"_",filesIn[,"name"],".",gr,".txt")
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
	write.xlsx2(stattab(httab,servRow),filexlsx,sheet=paste(sub("_mergedFeatures","",gr),"Identifed",sep="_"),append=TRUE)

	# NF<-max(colSums(httab))/colSums(httab)
	# norm2all<-httab*max(colSums(httab))/matrix(colSums(httab),byrow=TRUE,ncol=ncol(httab),nrow=nrow(httab))
	# c2<-colSums(httab[!(rownames(httab) %in% servRow),])
	# NF2<-max(c2)/c2
	# norm2feat<-httab*max(c2)/matrix(c2,byrow=TRUE,ncol=ncol(httab),nrow=nrow(httab))
	# write.xlsx2(stattab(norm2feat,servRow),filexlsx,sheet="Normalised to features mapped reads",append=TRUE)
	sampleTable<-cbind(Sample=colnames(httab),nr=unlist(GroupsSel[filesIn[,"rf"]]))
	rownames(sampleTable)<-colnames(httab)
	colData<-as.data.frame(sampleTable) # ,stringAsFactors=FALSE
	sets<-matrix(c("test","control"),nrow=2)
	for(s in 1:ncol(sets)){
		sel<-sampleTable[sampleTable[,"nr"] %in% sets[,s],"Sample"]
		out<-makeDG(httab,sel,colData,txt=paste(sub("_mergedFeatures","",gr),sep="_"),sets[1,s],sets[2,s],filexlsx)
	}

	if(nrow(htexp)>1) figVen(htexp,lim,paste(sub("_mergedFeatures","",gr),"venn diagram"),filexlsx)

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
			if((sum(is.na(c))==0 && length(table(c))>1) || ncol(c)<20){
				wb<-openxlsx::loadWorkbook(filexlsx)
				if(ncol(c)<20){
					png(paste0(fileName,"_corrplot.png"),width = 3+ncol(stat)/4, height = 1+ncol(stat)/4, res=150, units = "in")
					corrplot.mixed(c, p.mat = p, number.cex = 1, sig.level = .05,title=paste(sub("_mergedFeatures","",gr),set,method),mar=c(1,1,3,1))
					dev.off()
					openxlsx::insertImage(wb,sheet=substr(paste(sub("_mergedFeatures","",gr),set,method),0,31),file=paste0(fileName,"_corrplot.png"),width = 3+ncol(stat)/4, height = 1+ncol(stat)/4, dpi=150,startCol = 4)
					# openxlsx::insertPlot(wb,sheet=substr(paste(sub("_mergedFeatures","",gr),set,method),0,31),width = 3+ncol(stat)/4, height = 1+ncol(stat)/4, dpi=150,startCol = 4)
				}
				if(sum(is.na(c))==0 && length(table(c))>1){
					png(paste0(fileName,"_heatmap.png"),width = 3+ncol(stat)/4, height = 1+ncol(stat)/4, res=150, units = "in")
					heatmap.2(c,Rowv=TRUE,Colv=TRUE, dendrogram="row", revC = T, scale="none", col=greenred(75),na.rm=TRUE, key=TRUE, density.info="none", trace="none",mar=c(8,8))
					title(paste(sub("_mergedFeatures","",gr),set,method),cex.main=0.8)
					dev.off()
					openxlsx::insertImage(wb,sheet=substr(paste(sub("_mergedFeatures","",gr),set,method),0,31),file=paste0(fileName,"_heatmap.png"),width = 3+ncol(stat)/4, height = 1+ncol(stat)/4, dpi=150,startCol = 15)
					# openxlsx::insertPlot(wb,sheet=substr(paste(sub("_mergedFeatures","",gr),set,method),0,31),width = 3+ncol(stat)/4, height = 1+ncol(stat)/4, dpi=150,startCol = 15)
				}
				saveWorkbook(wb,filexlsx, overwrite = TRUE)
			}
		}
	}
	file.remove(paste0(fileName,"_corrplot.png"))
	file.remove(paste0(fileName,"_heatmap.png"))
}

warnings()
date()


