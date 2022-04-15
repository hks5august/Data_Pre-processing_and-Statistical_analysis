options(scipen=999)
library(MASS);
args <- commandArgs(TRUE)
library("caret")
train=read.csv(args[1],check.names=FALSE)
nzv <- nearZeroVar(train)
filteredDescr <- train[, -nzv]
#print dim(filteredDescr);
procValues <- preProcess(filteredDescr, method = c("center", "scale"))
scaledTraindata <-  predict(procValues, filteredDescr)
write.csv(scaledTraindata,file= "tt")
