#!/usr/bin/R --no-save
WD<-file.path(wd,"www","db","blast")
if(!dir.exists(WD)) dir.create(WD,recursive = TRUE, mode = "0777")
file.remove(file.path(WD,"db.done"))
system(paste("cd",WD," && update_blastdb --decompress nt"))
file.create(file.path(WD,"db.done"))
