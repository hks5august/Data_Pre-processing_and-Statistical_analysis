options(scipen=999)
library(MASS);
args <- commandArgs(TRUE)
library("caret")
mat=read.table(args[1],sep=",",check.names=FALSE, header=T, row.names=1)
nzv <- nearZeroVar(mat)
filteredDescr <- mat[, -nzv]
#print dim(filteredDescr);
write.table(filteredDescr,file= args[2],sep=",")
