setwd("C:/Users/Utente/Desktop/Andrea_ClonalEvolution/2 batch/parte A")
filelist<-list.files(pattern = ".txt")

reads<-lapply(filelist, read.delim, sep="\t", fill=TRUE, skip=581, nrows=2749694)
gc()

smoothlist<-lapply(reads,"[",  ,6)


names(smoothlist)<-filelist

for (i in 1:length(smoothlist)) {
  write.table(smoothlist[i], file=paste0(names(smoothlist)[i],"_SS.txt"), row.names=FALSE)
}
rm(list=ls())
setwd("C:/Users/Utente/Desktop/Andrea_ClonalEvolution/2 batch/parte B")

filelist<-list.files(pattern = ".txt")

reads<-lapply(filelist, read.delim, sep="\t", fill=TRUE, skip=691, nrows=2749694)
gc()

smoothlist<-lapply(reads,"[",  ,6)


names(smoothlist)<-filelist

for (i in 1:length(smoothlist)) {
  write.table(smoothlist[i], file=paste0(names(smoothlist)[i],"_SS.txt"), row.names=FALSE)
}
rm(list=ls())

setwd("C:/Users/Utente/Desktop/Andrea_ClonalEvolution/3 batch")
filelist<-list.files(pattern = ".txt")

reads<-lapply(filelist, read.delim, sep="\t", fill=TRUE, skip=691, nrows=2749694)
gc()

smoothlist<-lapply(reads,"[",  ,6)


names(smoothlist)<-filelist

for (i in 1:length(smoothlist)) {
  write.table(smoothlist[i], file=paste0(names(smoothlist)[i],"_SS.txt"), row.names=FALSE)
}
