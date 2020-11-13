#!/usr/local/bin/Rscript --vanilla
arg<-commandArgs()
WD<-file.path(arg[length(arg)],"www","db","blast")
if(!dir.exists(WD)) dir.create(WD,recursive = TRUE, mode = "0777")
file.remove(file.path(WD,"db.done"))
while(system(paste("cd",WD," && update_blastdb --decompress nt"))!=0) next;
file.create(file.path(WD,"db.done"))

