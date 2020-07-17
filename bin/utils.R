pkgs <- rownames(installed.packages())
file=paste0("~/R/pkgs_",R.Version()$major,".",R.Version()$minor,".txt")
write.table(pkgs[order(pkgs)], file=file,row.names=FALSE,col.names = FALSE,quote=FALSE)
Sys.setlocale(category = "LC_ALL", locale = "en_US.utf8")

Update<-function(x=""){
	if (!requireNamespace("BiocManager"))
		install.packages("BiocManager")
	chooseCRANmirror(graphics =FALSE,ind=1)
	chooseBioCmirror(graphics =FALSE,ind=1)
	BiocManager::install(character(), ask=FALSE)
	pkgs <- rownames(installed.packages())
	file=paste0("~/R/pkgs_",R.Version()$major,".",R.Version()$minor,".txt")
	write.table(pkgs[order(pkgs)], file=file,row.names=FALSE,col.names = FALSE,quote=FALSE)
	if(x!="") BiocManager::install(x)
	else {
		x<-select.list(c("Reinstall from list","Exit"),graphics =FALSE)
		if(x=="Reinstall from list"){
			sel<-select.list(dir("~/R",pattern = ".txt", recursive=FALSE),graphics =FALSE)
			pkgs<-read.table(file=paste0("~/R/",sel),as.is=TRUE)[,1]
			sp<-pkgs[!(pkgs %in% rownames(installed.packages()))]
			selp<-select.list(c("All_reinstall","All_packages","Exit",sp),graphics =FALSE)
			if(selp=="All_reinstall") BiocManager::install(pkgs, type="source")
			if(selp=="All_packages") BiocManager::install(sp, type="source")
			if(selp %in% pkgs) BiocManager::install(selp, type="source")
		}
	}
}

write.xlsx2<-function(data=c(),filexlsx="Book1.xlsx",sheet="Sheet1",append=FALSE,col.names=TRUE,row.names=TRUE){
	library(openxlsx)
	if(append & file.exists(filexlsx)){ wb <- loadWorkbook(filexlsx)
	} else wb<-createWorkbook()
	if(nchar(sheet)>31) sheet<-substr(sheet,1,31)
	addWorksheet(wb = wb, sheetName = sheet, gridLines = TRUE)
	freezePane(wb, sheet, firstRow = TRUE, firstCol = TRUE)
	writeData(wb, sheet = sheet, data, colNames = col.names, rowNames = row.names)
	saveWorkbook(wb,filexlsx,overwrite = TRUE)
}

RevCompl <- function(x) chartr("ATGCatgcUu","TACGtacgAa",sapply(lapply(strsplit(x,NULL),rev),paste,collapse=""))

