options(scipen=999)
library(MASS);
args <- commandArgs(trailingOnly = TRUE)
library(caret)
###### input file should be like samples in rows and genes should be in columns ####
train = read.csv(args[1],sep=",",header=T,row.names=1,check.names=FALSE)
test = read.csv(args[2],sep=",",header=T,row.names=1,check.names=FALSE)

train=data.frame(train)
test=data.frame(test)

nzv <- nearZeroVar(train)
filteredDescr <- train[, -nzv]
filteredDescr1 <- test[, -nzv]

pp <- preProcess(filteredDescr, method = c("range"))
tr_norm <-  predict(pp, filteredDescr)
te_norm <-  predict(pp, filteredDescr1)

#pp = preProcess(train, method = "range")
#pp
#tr_norm<-data.frame(predict(pp, train))
#tr_norm<-predict(pp, train)
#head(tr_norm)
tr_norm=data.frame(tr_norm)

#te_norm<-predict(pp, test)
#head(te_norm)
te_norm=data.frame(te_norm)

write.table(data.frame('Gene_ID'=rownames(tr_norm), tr_norm),file=args[3],sep=',',quote = F,row.names = FALSE)
write.table(data.frame('Gene_ID'=rownames(te_norm), te_norm),file=args[4],sep=',',quote = F,row.names = FALSE)
