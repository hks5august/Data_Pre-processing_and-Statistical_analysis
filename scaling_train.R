options(scipen=999)
library(MASS);
args <- commandArgs(TRUE)
library("caret")
train=read.csv(args[1],check.names=FALSE,row.names=1)
#test=read.csv(args[2],check.names=FALSE)
nzv <- nearZeroVar(train)
filteredDescr <- train[, -nzv]
#filteredDescr1 <- test[, -nzv]
#print dim(filteredDescr);
#scaledTestdata <-  predict(procValues, filteredDescr1)
write.csv(filteredDescr ,file= "tt")
#write.csv(scaledTestdata,file="ttt")
