#!/usr/bin/R --no-save
main_specie<-specie
DD<-file.path(wd,"www","db/")
archiveDir<-paste0(DD,"archive/")
main_specie_short<-paste0(substr(main_specie,1,1),substr(strsplit(main_specie,"_")[[1]][2],1,2))

# head -n 100000  $DV.gtf | grep "_biotype" | awk -F "_biotype" '{print $2}' | awk -F ";" '{print $1}' | sort | uniq -c | sort -n
faFile="piRBase_hsa.fa"; type="piRBase"
fa2gtf<-function(archiveDir,faFile="hg19-tRNAs.fa",outDir="gtf_biotypes/",type="GtRNAdb"){
	faTab<-paste0(archiveDir,faFile,".faTab")
	fasta2tab<-" | awk '/^>/ {printf(\"%s%s\\t\",(N>0?\"\\n\":\"\"),$1);N++;next;} {printf(\"%s\",$0);} END {printf(\"\\n\");}' > "
	system(paste0("cat ",archiveDir,faFile,fasta2tab,faTab))
	fa<-read.table(faTab,colClasses = "character",quote = "")
	for(d in fa[duplicated(fa[,ncol(fa)]),ncol(fa)]){
		dupl<-unique(fa[fa[,ncol(fa)]==d,1])
		if(length(dupl)>1)
			if(length(unique(substr(dupl,0,min(nchar(dupl))-1)))==1)
				fa[fa[,ncol(fa)]==d,1] <- paste0(unique(substr(dupl,0,min(nchar(dupl))-1)),paste(substr(dupl,min(nchar(dupl)),max(nchar(dupl))),collapse = "/"))
		else fa[fa[,ncol(fa)]==d,1] <- paste(unique(fa[fa[,ncol(fa)]==d,1]),collapse = "/")
	}
	fa<-fa[!duplicated(fa[,ncol(fa)]),]
	write.table(cbind(apply(cbind(">",fa[,1]),1,paste,collapse=""),fa[,ncol(fa)]),file=paste0(archiveDir,type,".fa"),sep="\n",col.names=FALSE,row.names=FALSE,quote=FALSE)
	system(paste0("cat ",archiveDir,type,".fa | bowtie -a --best -S --strata -p ",core," -v 0 -f ",DD,main_specie,"/",main_specie," - | ",samtools," view -F4 -o ",archiveDir,type,".sam") )
	sam<-read.table(paste0(archiveDir,type,".sam"),sep="\t",colClasses = "character")
	gtf <- cbind(sam[,3],type,"exon",sam[,4],as.integer(as.numeric(sam[,4])+nchar(sam[,10])),as.integer(sam[,5]),sam[,2],".",paste0("gene_id \"",sam[,1],"\"; transcript_id \"",sam[,1],"\";"))
	gtf[gtf[,7]=="0" ,7]<-"+"
	gtf[gtf[,7]=="16",7]<-"-"
	ch<-gtf[,1]
	ch[ch=="X"]<-"23"
	ch[ch=="Y"]<-"24"
	ch[ch=="MT"]<-"25"
	uc<-sort(unique(ch[!(ch %in% 1:25)]))
	for(c in 1:length(uc)) ch[ch==uc[c]] <- (c+25)
	gtf<-gtf[order(as.numeric(ch),as.numeric(gtf[,4])),]
	write.table(gtf,paste0(outDir,type,".gtf"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
}

gtfQuote<-function(file="gtf_biotypes/lncipedia_hc.gtf"){
	file.copy(file,sub(".gtf","_notQuoted.gtf",file),overwrite = TRUE)
	gtf<-read.table(file,sep="\t",colClasses = "character",quote = "")
	head(gtf)
	row<-1
	for(row in 1:nrow(gtf)){
		desc<-strsplit(gtf[row,9],";")
		type<-sub(" .*","",trimws(desc[[1]]))
		id<-sub("\"\"","\"",paste0("\"",sub(".* ","",trimws(desc[[1]])),"\""))
		gtf[row,9]<-paste(paste(type,id),collapse = "; ")
	}
	write.table(gtf,file=file,sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
}

gtfMergeFeatures<-function(file="gtf_biotypes/hg38_repeatmasker_2.gtf"){
	gtfAll<-read.table(file,sep="\t",colClasses = "character",quote = "")
	before<-nrow(gtfAll)
	gtfOut<-gtfAll[-c(1:before),]
	chs <- unique(gtfAll[,1])
	for(ch in chs){
		gtf<-gtfAll[gtfAll[,1]==ch,]
		row<-1
		while(row<=nrow(gtf)){
			overlap<-gtf[
				(as.numeric(gtf[,4])<=as.numeric(gtf[row,4]) & as.numeric(gtf[,5])>=as.numeric(gtf[row,4])) | 
				(as.numeric(gtf[,4])<=as.numeric(gtf[row,5]) & as.numeric(gtf[,5])>=as.numeric(gtf[row,5])) |
				(as.numeric(gtf[,4])>=as.numeric(gtf[row,4]) & as.numeric(gtf[,5])<=as.numeric(gtf[row,5])) ,]
			if(nrow(overlap)>1){
				desc<-strsplit(trimws(overlap[,9]),";")
				dTab<-rbind(sub(".* ","",trimws(desc[[1]])))
				colnames(dTab)<-sub(" .*","",trimws(desc[[1]]))
				chain<-""
				if(length(table(overlap[,7]))>1){
					gtf[row,7]<-"."
					chain<-paste0("[",overlap[,7],"]")
				}
				for(i in 2:length(desc)){
					dT2<-rbind(sub(".* ","",trimws(desc[[i]])))
					colnames(dT2)<-sub(" .*","",trimws(desc[[i]]))
					en<-colnames(dT2)[colnames(dT2) %in% colnames(dTab)]
					dTab<-rbind(dTab,NA)
					dTab[nrow(dTab),en]<-dT2[1,en]
					nn<-colnames(dT2)[!(colnames(dT2) %in% en)]
					dTab<-cbind(dTab,rbind(matrix(NA,nrow=nrow(dTab)-1,ncol=length(nn)),dT2[,nn]))
				}
				desc<-apply(dTab,2,function(x){ paste0("\"",paste(unique(paste0(gsub('"',"",x),chain)),collapse = "/"),"\""); }  )
				desc<- paste(paste0(paste(names(desc),desc),";"),collapse = " ")
				gtf[row,]<-c(as.character(gtf[row,1:3]),as.integer(min(as.numeric(overlap[,4]))),as.integer(max(as.numeric(overlap[,5]))),as.character(gtf[row,6:8]),desc)
				gtf<-gtf[!(rownames(gtf) %in% rownames(overlap[rownames(overlap)!=rownames(gtf)[row],])),]
			} else row <- row+1
		}
		gtfOut<-rbind(gtfOut,gtf)
	}
	write.table(gtfOut,file=sub(".gtf","_mergedFeatures.gtf",file),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
	after<-nrow(gtfOut)
	print(paste0(before,"-",after,"=",before-after," (",round((before-after)/before*100),"%) features merged to others"))
}

setwd(DD)
dir.create("gtf_biotypes")
dir.create(archiveDir)

download.file(paste0('http://regulatoryrna.org/database/piRNA/download/archive/v1.0/fasta/piR_human_v1.0.fa.gz'),paste0(archiveDir,"piRBase_hsa.fa.gz"),"auto",mode = "wb")
system(paste0("zcat ",archiveDir,"piRBase_hsa.fa.gz > ",archiveDir,"piRBase_hsa.fa") )
fa2gtf(archiveDir,"piRBase_hsa.fa","gtf_biotypes/","piRBase")
gtfMergeFeatures("gtf_biotypes/piRBase.gtf")

download.file(paste0('https://www.pirnadb.org/download/downloadarchive/gff_gtf/pirnadb.hg38.gtf.gz'),paste0(archiveDir,"pirnadb.hg38.gtf.gz"),"auto",mode = "wb")
download.file(paste0('https://www.pirnadb.org/download/downloadarchive/pirna/piRNAdb.hsa.v1_7_5.fa.gz'),paste0(archiveDir,"piRNAdb_hsa.gz"),"auto",mode = "wb")
system(paste0("zcat ",archiveDir,"piRNAdb_hsa.gz > ",archiveDir,"piRNAdb_hsa.fa") )
fa2gtf(archiveDir,"piRNAdb_hsa.fa","gtf_biotypes/","piRNAdb")
gtfMergeFeatures("gtf_biotypes/piRNAdb.gtf")

download.file(paste0('http://pirnabank.ibab.ac.in/downloads/all/human_all.zip'),paste0(archiveDir,"piRNAbank_hsa.zip"),"auto",mode = "wb")
system(paste0("zcat ",archiveDir,"piRNAbank_hsa.zip | sed '/^[^>]/ y/uU/tT/' | tr ':\n' '\t' | sed 's/Homo /homo_/g' | sed 's/\\t>/\\n>/g' > ",archiveDir,"piRNAbank_hsa.tsv"))
collapsed<-t(unique(read.table(paste0(archiveDir,"piRNAbank_hsa.tsv"))[,c(1,6)]))
cat(as.character(collapsed),file=paste0(archiveDir,"piRNAbank_hsa_collapsed.fa"),sep="\n")
fa2gtf(archiveDir,"piRNAbank_hsa_collapsed.fa","gtf_biotypes/","piRNAbank")
gtfMergeFeatures("gtf_biotypes/piRNAbank.gtf")

download.file(paste0('http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi19/hg19-tRNAs.fa'),paste0(archiveDir,"hg19-tRNAs.fa"),"auto",mode = "wb")
fa2gtf(archiveDir,"hg19-tRNAs.fa","gtf_biotypes/","GtRNAdb")
gtfMergeFeatures("gtf_biotypes/GtRNAdb.gtf")

system(paste0("tRNAscan-SE -qQ --detail -o# -m# -f# -l# -s# ",archiveDir,"GtRNAdb.fa"))
GtRNAdb<-read.table("gtf_biotypes/GtRNAdb.gtf",colClasses = "character",sep="\t",quote = "")
GtRNAfa<-read.table(paste0(archiveDir,"hg19-tRNAs.fa.faTab"),colClasses = "character",quote = "")
GtRNAfa<-GtRNAfa[!duplicated(GtRNAfa[,ncol(GtRNAfa)]),]
tRNAhalves<-tRF<-c()
i<-191
for(i in 1:nrow(GtRNAfa)){
	name<-GtRNAfa[i,1]
	str<-system(paste0("grep -A 5 -P \"",name,"[/.]\" ",archiveDir,"GtRNAdb.ss | grep ^Str | sed 's/^Str: //'"),intern=TRUE)
	len<-nchar(GtRNAfa[i,ncol(GtRNAfa)])
	tRFs<-as.numeric(gregexpr("..>\\.\\.\\.\\.\\.",str,fixed=FALSE)[[1]])+1
	# tRFs<-as.numeric(gregexpr(">>>...",str,fixed=TRUE)[[1]])+1
	tRFe<-as.numeric(gregexpr("\\.\\.\\.\\.\\.<..",str,fixed=FALSE)[[1]])+6
	if(length(tRFs)<3 & length(tRFe)<3){
		tRFs<-as.numeric(gregexpr("..>\\.\\.\\.\\.",str,fixed=FALSE)[[1]])[c(1,2,4)]+1
		tRFe<-as.numeric(gregexpr("\\.\\.\\.\\.<..",str,fixed=FALSE)[[1]])[c(1,2,4)]+5
	}
	# if(length(tRFs)>3 & length(tRFe)>3){
	# 	tmp<-!((tRFe-tRFs)==min(tRFe-tRFs))
	# 	tRFe<-tRFe[tmp]
	# 	tRFs<-tRFs[tmp]
	# }
	if(length(tRFs)!=3 | length(tRFe)!=3){ print(i);  }
	# substr(GtRNAfa[i,13],tRFs[2],tRFe[2])
	j<-1
	for(j in grep(paste0(name,"\""),GtRNAdb[,9])){
		tRFpre5<-tRF5D<-tRFi<-tRF3T<-tRF1<-half3<-half5<-tmp<-GtRNAdb[j,]
		tRFpre5[9] <- gsub(name,paste0(name,"-tRFpre5"),tmp[9])
		tRF5D[9]   <- gsub(name,paste0(name,"-tRF5D"  ),tmp[9])
		tRFi[9]    <- gsub(name,paste0(name,"-tRFi"   ),tmp[9])
		tRF3T[9]   <- gsub(name,paste0(name,"-tRF3T"  ),tmp[9])
		tRF1[9]    <- gsub(name,paste0(name,"-tRF1"   ),tmp[9])
		half3[9]   <- gsub(name,paste0(name,"-half3"  ),tmp[9])
		half5[9]   <- gsub(name,paste0(name,"-half5"  ),tmp[9])
		if(tmp[7]=="+"){
			tRFpre5[4]<- as.numeric(tmp[4])-20
			tRFpre5[5]<- as.numeric(tmp[4])-4
			tRF5D[5]  <- as.numeric(tmp[4])+tRFs[1]
			tRFi[4]   <- as.numeric(tmp[4])+tRFe[1]
			tRFi[5]   <- as.numeric(tmp[4])+tRFs[3]
			tRF3T[4]  <- as.numeric(tmp[4])+tRFe[3]
			tRF1[4]   <- as.numeric(tmp[5])+4
			tRF1[5]   <- as.numeric(tmp[5])+20
			half3[4]  <- as.numeric(tmp[4])+tRFe[2]
			half5[5]  <- as.numeric(tmp[4])+tRFs[2]
		} else {
			tRFpre5[4]<- as.numeric(tmp[5])+4
			tRFpre5[5]<- as.numeric(tmp[5])+20
			tRF5D[4]  <- as.numeric(tmp[5])-tRFs[1]
			tRFi[4]   <- as.numeric(tmp[5])-tRFs[3]
			tRFi[5]   <- as.numeric(tmp[5])-tRFe[1]
			tRF3T[5]  <- as.numeric(tmp[5])-tRFe[3]
			tRF1[4]   <- as.numeric(tmp[4])-20
			tRF1[5]   <- as.numeric(tmp[4])-4
			half3[5]  <- as.numeric(tmp[5])-tRFe[2]
			half5[4]  <- as.numeric(tmp[5])-tRFs[2]
		}
		tRF<-rbind(tRF,tRFpre5,tRF5D,tRFi,tRF3T,tRF1)
		tRNAhalves<-rbind(tRNAhalves,half3,half5)
	}
}
ch<-tRF[,1]
ch[ch=="X"]<-"23"
ch[ch=="Y"]<-"24"
tRF<-tRF[order(as.numeric(ch),as.numeric(tRF[,4])),]
write.table(tRF,"gtf_biotypes/tRF.gtf",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
gtfMergeFeatures("gtf_biotypes/tRF.gtf")
# system(paste0("$HOME/data/conda/bin/tabix -p gff gtf_biotypes/tRF.gtf"),intern=TRUE)
ch<-tRNAhalves[,1]
ch[ch=="X"]<-"23"
ch[ch=="Y"]<-"24"
tRNAhalves<-tRNAhalves[order(as.numeric(ch),as.numeric(tRNAhalves[,4])),]
write.table(tRNAhalves,"gtf_biotypes/tRNAhalves.gtf",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
gtfMergeFeatures("gtf_biotypes/tRNAhalves.gtf")

download.file(paste0('ftp://mirbase.org/pub/mirbase/CURRENT/genomes/',main_specie_short,".gff3"),paste0(archiveDir,"miRBase_",main_specie_short,".gff"),"auto",mode = "wb")
library("rtracklayer")
tmp<-import.gff3(paste0(archiveDir,"miRBase_",main_specie_short,".gff"))
export.gff2(tmp,paste0(archiveDir,"miRBase_",main_specie_short,".gtf"))
system(paste0("cat ",archiveDir,"miRBase_",main_specie_short,".gtf | grep '\tmiRNA\t' | sed 's/^chr//' | sed 's/\tmiRNA\t/\texon\t/' | sed 's/; Name /; gene_id /' | sed 's/\tID /\ttranscript_id /' | sed 's/$/;/' | sed 's/;;/;/' > gtf_biotypes/miRBase_mature.gtf"),intern=TRUE)
gtfMergeFeatures("gtf_biotypes/miRBase_mature.gtf")
system(paste0("cat ",archiveDir,"miRBase_",main_specie_short,".gtf | grep '\tmiRNA_primary_transcript\t' | sed 's/^chr//' | sed 's/\tmiRNA_primary_transcript\t/\texon\t/' | sed 's/; Name /; gene_id /' | sed 's/\tID /\ttranscript_id /' > gtf_biotypes/miRBase_hairpin.gtf"),intern=TRUE)
gtfMergeFeatures("gtf_biotypes/miRBase_hairpin.gtf")
		
download.file(paste0('https://lncipedia.org/downloads/lncipedia_5_2/full-database/lncipedia_5_2_hg38.gtf'),paste0(archiveDir,"lncipedia_",main_specie_short,".gtf"),"auto",mode = "wb")
system(paste0("cat ",archiveDir,"lncipedia_",main_specie_short,".gtf | grep -v '^track' | sed 's/^chr//' > gtf_biotypes/lncipedia.gtf"))
gtfQuote("gtf_biotypes/lncipedia.gtf")
gtfMergeFeatures("gtf_biotypes/lncipedia.gtf")
download.file(paste0('https://lncipedia.org/downloads/lncipedia_5_2/high-confidence-set/lncipedia_5_2_hc_hg38.gtf'),paste0(archiveDir,"lncipedia_hc_",main_specie_short,".gtf"),"auto",mode = "wb")
system(paste0("cat ",archiveDir,"lncipedia_hc_",main_specie_short,".gtf | grep -v '^track' | sed 's/^chr//' > gtf_biotypes/lncipedia_hc.gtf"))
gtfQuote("gtf_biotypes/lncipedia_hc.gtf")
gtfMergeFeatures("gtf_biotypes/lncipedia_hc.gtf")

system(paste0("cat ",DD,main_specie,"/",main_specie,".gtf | grep _biotype | awk -F \"_biotype\" '{print $2}' | awk -F \";\" '{print $1}' | sort | uniq -c | sort -rn > gtf_biotypes/",main_specie,".gtf_biotypes.txt"))
biotypes<-c("rRNA","snoRNA","miRNA","snRNA","misc_RNA","Mt_tRNA","Mt_rRNA","protein_coding","lncRNA","processed_pseudogene","vaultRNA")
for(type in biotypes){
	system(paste0("grep '_biotype \"",type,"\"' ",DD,main_specie,"/",main_specie,".gtf | grep '\texon\t' > gtf_biotypes/",type,".gtf"))
	gtfMergeFeatures(paste0("gtf_biotypes/",type,".gtf"))
}
system(paste0("grep ^MT ",DD,main_specie,"/",main_specie,".gtf | grep '\texon\t' > gtf_biotypes/MT.gtf"))
gtfMergeFeatures("gtf_biotypes/MT.gtf")
system(paste0("grep -v \"Mt_tRNA\\|Mt_rRNA\" gtf_biotypes/MT.gtf | grep '\texon\t' > \"gtf_biotypes/other(MT).gtf\""))
gtfMergeFeatures("gtf_biotypes/other(MT).gtf")
system(paste0("grep \"Y_RNA\\|RNY\" gtf_biotypes/misc_RNA.gtf | grep '\texon\t' > \"gtf_biotypes/YRNA(misc_RNA).gtf\""))
gtfMergeFeatures("gtf_biotypes/YRNA(misc_RNA).gtf")
system(paste0("grep pRNA gtf_biotypes/misc_RNA.gtf | grep '\texon\t' > \"gtf_biotypes/pRNA(misc_RNA).gtf\""))
gtfMergeFeatures("gtf_biotypes/pRNA(misc_RNA).gtf")
system(paste0("grep -v \"Y_RNA\\|RNY\\|pRNA\" gtf_biotypes/misc_RNA.gtf | grep '\texon\t' > \"gtf_biotypes/noYorPiwi(misc_RNA).gtf\""))
gtfMergeFeatures("gtf_biotypes/noYorPiwi(misc_RNA).gtf")
system(paste0("grep '\tgene\t' ",DD,main_specie,"/",main_specie,".gtf | sed 's/\tgene_id \"/\tgene_id \"intron_/' | sed 's/\tgene\t/\texon\t/' > gtf_biotypes/Ensembl_genes.gtf"))
gtfMergeFeatures("gtf_biotypes/Ensembl_genes.gtf")

file.create("gtf_biotypes/All_catalogue_types.gtf")
for(type in c(biotypes,"MT")){
	system(paste0("cat gtf_biotypes/",type,".gtf >> gtf_biotypes/All_catalogue_types.gtf"))
}
system(paste0("cat gtf_biotypes/All_catalogue_types.gtf | sort | uniq > gtf_biotypes/All_catalogue_types_sorted.gtf"))
system(paste0("grep '\texon\t' ",DD,main_specie,"/",main_specie,".gtf | grep -v -xFf gtf_biotypes/All_catalogue_types_sorted.gtf > gtf_biotypes/Other_types.gtf"))
gtfMergeFeatures("gtf_biotypes/Other_types.gtf")




download.file("https://genome.ucsc.edu/cgi-bin/hgTables?clade=mammal&org=Human&db=hg38&hgta_group=rep&hgta_track=rmsk&hgta_table=rmsk&hgta_regionType=genome&hgta_outputType=gff&hgta_outFileName=hg38_RepeatMasker_UCSC.gtf.gz&hgta_compressType=gzip&hgta_doTopSubmit=get",paste0(archiveDir,"hg38_RepeatMasker_UCSC.gtf.gz"),"auto",mode = "wb")
system(paste0("zcat ",archiveDir,"hg38_RepeatMasker_UCSC.gtf.gz | sed 's/^chr//' > gtf_biotypes/hg38_repeatmasker.gtf"))
type<-"rRNA"
for(type in c("rRNA","tRNA")){
	system(paste0("grep '",type,"' gtf_biotypes/hg38_repeatmasker.gtf | grep '\texon\t' > gtf_biotypes/RepeatMasker_",type,".gtf"))
	gtfMergeFeatures(paste0("gtf_biotypes/RepeatMasker_",type,".gtf"))
}

unlink("gtf_biotypes/repeatMasker",recursive = TRUE)
dir.create("gtf_biotypes/repeatMasker")
chRs<-system("awk '{print $1}' gtf_biotypes/hg38_repeatmasker.gtf | sort | uniq",intern = TRUE)

library(foreach)
library(doMC)
registerDoMC()

# for(chR in chRs[order(chRs,decreasing = TRUE)]){
err<-foreach(chR=chRs[order(chRs,decreasing = TRUE)]) %dopar% {
	if(!file.exists(paste0("gtf_biotypes/repeatMasker/",chR,".done"))){
		system(paste0("grep '^",chR,"\t' gtf_biotypes/hg38_repeatmasker.gtf > gtf_biotypes/repeatMasker/",chR,".gtf"),intern = TRUE)
		gtfMergeFeatures(paste0("gtf_biotypes/repeatMasker/",chR,".gtf"))
		file.create(paste0("gtf_biotypes/repeatMasker/",chR,".done"))
	}
}
done<-sub(".done","",dir("gtf_biotypes/repeatMasker",pattern = ".done"))
chRs[!(chRs %in% done)]
file.create(paste0("gtf_biotypes/RepeatMasker_mergedFeatures.gtf"))
for(chR in chRs){
	if(file.exists(paste0("gtf_biotypes/repeatMasker/",chR,".done"))){
		system(paste0("grep '^",chR,"\t' gtf_biotypes/repeatMasker/",chR,"_mergedFeatures.gtf >> gtf_biotypes/RepeatMasker_mergedFeatures.gtf"),intern = TRUE)
	}
}
