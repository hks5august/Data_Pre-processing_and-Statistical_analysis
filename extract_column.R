filesT <- list.files(pattern="T.txt")
filesN <- list.files(pattern="N.txt")


genes <- read.table(filesT[1], header=T, sep="\t")[,1]     # gene names
dfT    <- do.call(cbind,lapply(filesT,function(fn)read.table(fn,header=T, sep="\t")[,2]))
dfN    <- do.call(cbind,lapply(filesN,function(fn)read.table(fn,header=T, sep="\t")[,2]))

genes<-as.data.frame(genes)
dfT<-as.data.frame(dfT)
dfN<-as.data.frame(dfN)

dfTN  <- cbind(genes,dfT,dfN)

write.table(dfTN,file="dfTN", quote = FALSE, sep = "\t")
