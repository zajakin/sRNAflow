library(gplots)
library(corrplot)
# library(edgeR)
library(DESeq2)
library(EnhancedVolcano) # sudo apt install libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libproj-dev
library(PCAtools)
library(gridExtra)
library(ggplot2)
library(VennDiagram)
library(openxlsx)
tmp<-futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

write2xlsx<-function(data=c(),wb,sheet="Sheet1",col.names=TRUE,row.names=TRUE){
	if(nchar(sheet)>31) sheet<-substr(sheet,1,31)
	addWorksheet(wb = wb, sheetName = sheet, gridLines = TRUE)
	freezePane(wb, sheet, firstRow = TRUE, firstCol = TRUE)
	writeData(wb, sheet = sheet, data, colNames = col.names, rowNames = row.names)
}

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
	wb<-createWorkbook()
	filexlsx<- file.path(ED,"isomiR-SEA",paste0(sample,"_isomiR-SEA.xlsx"))
	tab<-read.table(file.path(p,"summary.txt"), sep = " ",header = F)[,-2]
	write2xlsx(tab,wb,sheet="Summary",row.names = FALSE,col.names = FALSE)
	for(t in tabs){
		infile<-file.path(p,paste0("out_result_mature_",t,".txt"))
		if(length(readLines(infile,n=2))==2){
			tab<-read.table(infile, sep = "\t",header = F,skip = 1,comment.char = "",fill = T)
			colnames(tab)<-read.table(infile,header = FALSE,nrows = 1,skip = 0,comment.char = "")
			write2xlsx(tab,wb,sheet=t,row.names = FALSE)
		}
	}
	saveWorkbook(wb,filexlsx,overwrite = TRUE)
}
zip(file.path(ED,paste0(Exp,"_isomiR-SEA.zip")),files=dir(file.path(ED,"isomiR-SEA"),"_isomiR-SEA.xlsx",full.names = T,recursive = TRUE),flags="-joq9")

figVen<-function(dat,lim,limS,txt="",wbF){
	title<-paste0(txt," >",lim," in ",limS," samples")
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
		write2xlsx(delist,wbF,sheet=txt,row.names = F)
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
		insertPlot(wbF,sheet=txt,width = 8, height = 6, dpi=150,startCol = 4)
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

# myEdgeR<-function(qlf,sheet,y,wb){
# 	res<-topTags(qlf, n= nrow(y$counts) )
# 	# res<-res[as.numeric(res[[1]][,"PValue"])<1,]
# 	tmp<-(y$counts[rownames(res),]*rep(as.numeric(y$samples[,"norm.factors"]),each=nrow(res)))
# 	# colSums(tmp)
# 	# y$samples
# 	res<-cbind(rownames(res),res,NA,baseMedian=rowMedians(tmp),baseMean=rowMeans(tmp),NA,
# 		tmp,NA,y$counts[rownames(res),])
# 	write2xlsx(as.data.frame(res),wb,row.names=TRUE,col.names = TRUE,sheet=sheet)
# }
 
makeDG<-function(httab,sel,colData,txt="test",cat1=sets[1,s],cat2=sets[2,s],wb,servRow=c()){
	dat<-httab[!(rownames(httab) %in% servRow),sel]
	dat<-as.matrix(dat[rowSums(dat>lim)>limS,])
	mode(dat)<-"integer"
	# group<-as.character(colData[sel,"nr"])
	# design <- model.matrix(~group)
	colData[,2]<-as.factor(colData[,2])
	if(sum((rowSums((dat==0)+0)==0)+0)>1){
		dds <- DESeqDataSetFromMatrix(countData = dat, colData = colData[sel,], design = ~ nr)
		err<-try(dds2 <- DESeq(dds,quiet = T))
		if(typeof(err)!="S4"){
			write2xlsx(gsub("\n"," ",as.character(err)),wb,row.names=TRUE,col.names = TRUE,sheet=paste(txt,"DESeq2"))
			return(rbind(list(baseMean=1,log2FoldChange=1,lfcSE=1,stat=1,pvalue=1,padj=1))[-1,])
		}
		res <- as.data.frame(results(dds2,contrast=c("nr","test","control"),lfcThreshold=log2FoldChange, alpha=padj, cooksCutoff=TRUE))
		res<-res[order(res[,"pvalue"]),]
		out<-res
		res<-res[!is.na(res[,"padj"]) & res[,"padj"]<padj,]
		if(nrow(res)>0){
			counts<-rbind(as.matrix(counts(dds2,normalized=TRUE))[rownames(res),])
			res<-cbind(rbind(res)," "=" ",baseMedian=rowMedians(counts),"  "=" ",counts,"   "=" ",Raw=rbind(dat[rownames(res),]))
			write2xlsx(as.data.frame(rbind(res,NA)),wb,row.names=TRUE,col.names = TRUE,sheet=paste(txt,"DESeq2"))
			#TODO Make working for mouse etc...  ####
			# if(txt == "all_miRNA") try(write2xlsx(tomieaa(sub("\\[[+-.]\\]","",unlist(strsplit(sub(".*_mergedFeatures_","",rownames(res)),"/"))),analysis="ora"),
			# 									   wb,row.names=TRUE,col.names = TRUE,sheet=paste(txt,"MIEAA_ORA")))
			# if(txt == "all_miRNA") try(write2xlsx(tomieaa(sub("\\[[+-.]\\]","",unlist(strsplit(sub(".*_mergedFeatures_","",rownames(res)),"/"))),analysis="gsea"),
			# 									   wb,row.names=TRUE,col.names = TRUE,sheet=paste(txt,"MIEAA_GSEA")))
			#TODO Diff expression DESeq2  GOstats ####
			# if(txt %in% c("Ensembl_genes","protein_coding","all_lncRNA","all_miRNA")){
			# # 	#  https://www.bioconductor.org/packages/release/bioc/manuals/GOstats/man/GOstats.pdf
			# 	myGO(myids=sub("\\[[+-.]\\]","",unlist(strsplit(sub(".*_mergedFeatures_","",rownames(res)),"/"))), minimum.population=5,deUp=c(),deDown=c())
			# }
		} else write2xlsx(c("No matching results found"),wb,row.names=TRUE,col.names = TRUE,sheet=paste(txt,"DESeq2"))
	} else write2xlsx(c("Normalisation not possible: No found features presented in all samples"),wb,row.names=TRUE,col.names = TRUE,sheet=paste(txt,"DESeq2"))
	subtitle<-paste(txt,"significant = ",nrow(res))
	# writeData(wb2, sheet = txt, subtitle,startCol=1,startRow=2)
	# writeData(wb2, sheet = txt, t(colData),startCol=ncol(res)+4,colNames = F, rowNames = T,startRow=1)
	if(nrow(out)<500){
		print(EnhancedVolcano(out[,1:6],sub("^hsa-","",sub(".*_mergedFeatures_","",rownames(out))),"log2FoldChange","padj",pCutoff=padj, FCcutoff=log2FoldChange,
							  xlim = c(min(out[,"log2FoldChange"], na.rm = TRUE) - 0.5, max(out[,"log2FoldChange"], na.rm = TRUE) + 0.5),
							  ylim = c(0, max(-log10(out[,"padj"]), na.rm = TRUE) + 0.5),
							  col = c("grey30", "blue", "forestgreen", "red"),colAlpha=1,legendLabels = c("NS", expression(abs(Log[2] ~ FC) ), expression("adj. p-value"), expression(abs(Log[2] ~ FC) ~ " and adj. p-value")),
							  drawConnectors = T,maxoverlapsConnectors = 42,lengthConnectors=unit(0.01, "npc"),legendPosition = "top",pointSize = 1,subtitle=subtitle ))
		insertPlot(wb,sheet=paste(txt,"DESeq2"),width = 12, height = 12, dpi=150,startCol = 20,startRow = 5)
	}
	print(EnhancedVolcano(out[,1:6],NA,"log2FoldChange","padj",pCutoff=padj, FCcutoff=log2FoldChange,
						  xlim = c(min(out[,"log2FoldChange"], na.rm = TRUE) - 0.5, max(out[,"log2FoldChange"], na.rm = TRUE) + 0.5),
						  ylim = c(0, max(-log10(out[,"padj"]), na.rm = TRUE) + 0.5),
						  col = c("grey30", "blue", "forestgreen", "red"),colAlpha=1,legendLabels = c("NS", expression(abs(Log[2] ~ FC) ), expression("adj. p-value"), expression(abs(Log[2] ~ FC) ~ " and adj. p-value")),
						  drawConnectors = F,lengthConnectors=unit(0.01, "npc"),legendPosition = "top",pointSize = 1,subtitle=subtitle))
	insertPlot(wb,sheet=paste(txt,"DESeq2"),width = 9, height = 7, dpi=300,startCol = 10,startRow = 5)
	if(nrow(dat)>1){
		p <- as.matrix(counts(dds2,normalized=TRUE))
		rownames(p)<-sub("^hsa-","",sub(".*_mergedFeatures_","",rownames(p)))
		p <- pca(p, metadata = colData[sel,]) ## , removeVar = 0.1 -- removing the lower 10% of variables based on variance
		print(screeplot(p, axisLabSize = 18, titleLabSize = 22))
		insertPlot(wb,sheet=paste(txt,"DESeq2"),width = 9, height = 7, dpi=150,startCol = 10,startRow = 75)
		print(biplot(p, colby = 'nr', colLegendTitle = txt, encircle = TRUE, encircleFill = TRUE, hline = 0, vline = c(-25, 0, 25),
					 legendPosition = 'top', legendLabSize = 16, legendIconSize = 8.0, showLoadings = TRUE, sizeLoadingsNames = 5))
		insertPlot(wb,sheet=paste(txt,"DESeq2"),width = 9, height = 7, dpi=300,startCol = 10,startRow = 40)
	}
	# indat <- DGEList( counts = dat, group = group)
	# RLE <- calcNormFactors(indat,method="RLE")$samples[,"norm.factors"]
	# RLEgenes<-table(rowSums(dat>0)==ncol(dat))["TRUE"]
	# y <- calcNormFactors(indat,method="TMM")
	# corLibSize<-round(colSums(y$counts*rep(as.numeric(y$samples[,"norm.factors"]),each=nrow(y$counts))))
	# NFcorLibSize<-colSums(y$counts*rep(NF[sel],each=nrow(y$counts)))
	# RLEcorLibSize<-colSums(y$counts*rep(RLE,each=nrow(y$counts)))
	# foundGenes<-colSums(y$counts!=0)
	# write2xlsx(t(cbind(y$samples,"TMM corrected libSize"=corLibSize,
	# 		"Number of genes used for RLE(DESeq2) normalisation"=RLEgenes,
	# 		"RLE NF"=RLE,"NF by mapped"=NF[sel],"NF2 by mapped to features"=NF2[sel],
	# 		"LibSize corrected by RLE"=RLEcorLibSize,"LibSize corrected by NF"=NFcorLibSize,
	# 		foundGenes=foundGenes)),wb,sheet=paste(txt,"NF"))
	# y <- estimateDisp(y,design)
	# myEdgeR(glmQLFTest(glmQLFit(y,design),coef=2),paste(txt,"edgeR TMM quasi-likelihood F-tests"),y,wb)
	# myEdgeR(glmLRT(glmFit(y,design),coef=2),paste(txt,"edgeR TMM likelihood ratio tests"),y,wb)
	# indat$samples[,"norm.factors"]<-NF[sel]
	# y<-estimateDisp(indat,design)
	# myEdgeR(glmQLFTest(glmQLFit(y,design),coef=2),paste(txt,"edgeR LibSize quasi-likelihood F-tests"),y,wb)
	# indat$samples[,"norm.factors"]<-NF2[sel]
	# y<-estimateDisp(indat,design)
	# myEdgeR(glmQLFTest(glmQLFit(y,design),coef=2),paste(txt,"mappedFeat edgeR quasi-likelihood F-tests"),y,wb)
	# myEdgeR(glmLRT(glmFit(y,design),coef=2),paste(txt,"edgeR LibSize likelihood ratio tests"),y,wb)
	# return(out)
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
wb<-createWorkbook()
value<-c(Exp=Exp,specie=specie,tsize=tsize,Rep=Rep,blast=blast,ad3=ad3,ad5=ad5,sizerange =sizerange,lim=lim,limS=limS,log2FoldChange=log2FoldChange,padj =padj,email =email,smtpServer=smtpServer)
value<-rbind(cbind(variable=names(value),value=value," "=" ","  "=" ")," ",c("file","size","date","group"),cbind(FilesIn,group=GroupsSel[FilesIn[,"file"]]))
write2xlsx(value,wb,sheet="Settings",row.names = FALSE)

tax<-""
if(strategy== "successively") tax<-"_tax_filtRC"

stats<-unlist(lapply(dir(ED,full.names =TRUE),dir,pattern=paste0("stat",tax,".txt"),full.names =TRUE))
#system(paste0("paste ",paste(stats,collapse = " ")," | tr ' ' '\t' | cut -f ",paste(c(1,seq(2,length(stats)*2+1,2)),collapse = ","), "> ",ED,"/stats",tax,"__.tsv"),intern = TRUE)
rows <- grep("^Ensembl_genes_mergedFeatures",readLines(stats[1]))-1
stat<-data.frame(read.table(stats[1],na.strings = "",header = TRUE,nrows = rows,check.names = FALSE),row.names = 1,check.names = FALSE)
for(i in stats[-1]) stat<-cbind(stat,data.frame(read.table(i,na.strings = "",header = TRUE,nrows = rows,check.names = FALSE),row.names = 1,check.names = FALSE))
write.table(stat,paste0(ED,"/stats",tax,".tsv"),quote=FALSE, sep="\t",col.names=TRUE)
#stat<-read.table(file.path(ED,paste0("stats",tax,".tsv")), sep = "\t",header = TRUE)
enrow<-grep("all_Ensembl",rownames(stat))
write2xlsx(stat[1:enrow,],wb,sheet="Quality")
write2xlsx(stat[(enrow+1):nrow(stat),],wb,sheet="Catalog")

endrow<- grep("Ensembl_genes_mergedFeatures",rownames(stat))-1
samples <- rep(colnames(stat),each=(endrow-enrow))
RNA_types <- rep(rownames(stat)[(enrow+1):endrow], ncol(stat))
frequency <- as.numeric(unlist(stat[(enrow+1):endrow,]))
data <- data.frame(samples,RNA_types,frequency)
# png(paste0(fileName,".png"),width = 6+ncol(stat)/4, height = 8, res=300, units = "in")
print(ggplot(data, aes(fill=RNA_types, y=frequency, x=samples)) + geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
# dev.off()
insertPlot(wb,sheet="Catalog",width = 6+ncol(stat)/4, height = 8, dpi=300,startCol = 8)
# insertImage(wb,sheet="Catalog",file=paste0(fileName,".png"),width = 6+ncol(stat)/4, height = 8, dpi=300,startCol = 8)

if(exists("species99")){
	httab<-data.frame(species99,row.names = 1)
	for(i in 1:nrow(filesIn)){
		tmp<-data.frame(read.table(paste0(filesIn[i,"wd"],"forKrona/",filesIn[i,"name"],".counts",tax,".txt")),row.names = 2)
		httab[[filesIn[i,"name"]]]<- tmp[rownames(httab),]
	}
	httab[is.na(httab)]<-0
	httab<-httab[rowSums(httab[,-1])>0,]
	httab<-httab[order(-rowSums(httab[,-1])),]
	write2xlsx(httab,wb,sheet="Species")
}
saveWorkbook(wb,filexlsx, overwrite = TRUE)
# file.remove(paste0(fileName,".png"))
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

	htexp<-httab[rowSums(httab>lim)>limS,]
	htexp<-rbind(htexp[!(rownames(htexp) %in% servRow),])
	print(paste(gr,grep(gr,deGTF),"/",length(deGTF),"rows =",nrow(htexp)))
	if(nrow(htexp)==0) next;
	write2xlsx(stattab(httab,servRow),wb,sheet=paste(sub("_mergedFeatures","",gr),"Identifed",sep="_"))

	# NF<-max(colSums(httab))/colSums(httab)
	# norm2all<-httab*max(colSums(httab))/matrix(colSums(httab),byrow=TRUE,ncol=ncol(httab),nrow=nrow(httab))
	# c2<-colSums(httab[!(rownames(httab) %in% servRow),])
	# NF2<-max(c2)/c2
	# norm2feat<-httab*max(c2)/matrix(c2,byrow=TRUE,ncol=ncol(httab),nrow=nrow(httab))
	# write2xlsx(stattab(norm2feat,servRow),wb,sheet="Normalised to features mapped reads")
	sampleTable<-cbind(Sample=colnames(httab),nr=unlist(GroupsSel[filesIn[,"rf"]]))
	rownames(sampleTable)<-colnames(httab)
	colData<-as.data.frame(sampleTable) # ,stringAsFactors=FALSE
	sets<-matrix(c("test","control"),nrow=2)
	for(s in 1:ncol(sets)){
		sel<-sampleTable[sampleTable[,"nr"] %in% sets[,s],"Sample"]
		makeDG(httab,sel,colData,txt=paste(sub("_mergedFeatures","",gr),sep="_"),sets[1,s],sets[2,s],wb)
	}

	if(nrow(htexp)>1) figVen(htexp,lim,limS,paste(sub("_mergedFeatures","",gr),"venn diagram"),wb)

	for(set in c("expressed only")){ #"all", 
		if(set=="all"){ d<-httab[!(rownames(httab) %in% servRow),] } else d<-htexp
		d<-d[,colSums(d)>0]
		if(nrow(d)<3) next
		for(method in c("spearman")){ # "pearson", 
			c<-cor(d,use="pairwise.complete.obs",method=method)
			p <- cor.mtest(d, conf.level = .95)$p
			pv<-data.frame(signif(p,2),row.names = paste0(rownames(c),"_p"),check.names = F)
			colnames(pv)<-colnames(c)
			sheet<-substr(paste(sub("_mergedFeatures","",gr),set,method),0,31)
			write2xlsx(rbind(data.frame(signif(c,2),check.names = F)," "=" ","p-values"=colnames(c),pv),wb,sheet=sheet)
			if((sum(is.na(c))==0 && length(table(c))>1) || ncol(c)<20){
				if(ncol(c)<20){
					# png(paste0(fileName,"_corrplot.png"),width = 3+ncol(stat)/4, height = 1+ncol(stat)/4, res=150, units = "in")
					corrplot.mixed(c, p.mat = p, number.cex = 1, sig.level = .05,title=paste(sub("_mergedFeatures","",gr),set,method),mar=c(1,1,3,1))
					# dev.off()
					# insertImage(wb,sheet=sheet,file=paste0(fileName,"_corrplot.png"),width = 3+ncol(stat)/4, height = 1+ncol(stat)/4, dpi=150,startCol = 4)
					insertPlot(wb,sheet=sheet,width = 3+ncol(stat)/4, height = 3+ncol(stat)/4, dpi=150,startCol = 4)
				}
				if(sum(is.na(c))==0 && length(table(c))>1){
					# png(paste0(fileName,"_heatmap.png"),width = 3+ncol(stat)/4, height = 3+ncol(stat)/4, res=150, units = "in")
					heatmap.2(c,Rowv=TRUE,Colv=TRUE, dendrogram="row", revC = T, scale="none", col=greenred(75),na.rm=TRUE, key=TRUE, density.info="none", trace="none",mar=c(8,8))
					title(paste(sub("_mergedFeatures","",gr),set,method),cex.main=0.8)
					# dev.off()
					# insertImage(wb,sheet=sheet,file=paste0(fileName,"_heatmap.png"),width = 3+ncol(stat)/4, height = 3+ncol(stat)/4, dpi=150,startCol = 15)
					insertPlot(wb,sheet=sheet,width = 3+ncol(stat)/4, height = 3+ncol(stat)/4, dpi=150,startCol = 15)
				}
				# saveWorkbook(wb,filexlsx, overwrite = TRUE)
			}
		}
	}
	# file.remove(paste0(fileName,"_corrplot.png"))
	# file.remove(paste0(fileName,"_heatmap.png"))
}
saveWorkbook(wb,filexlsx, overwrite = TRUE)

warnings()
date()


