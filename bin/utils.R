Sys.setlocale(category = "LC_ALL", locale = "en_US.utf8")

Update<-function(x=""){
	if (!requireNamespace("BiocManager"))
		install.packages("BiocManager")
	chooseCRANmirror(graphics =FALSE,ind=1)
	chooseBioCmirror(graphics =FALSE,ind=1)
	BiocManager::install(character(), ask=FALSE)
	if(x!="") BiocManager::install(x)
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

