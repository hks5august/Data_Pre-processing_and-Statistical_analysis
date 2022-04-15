####Usage Rscript quantile_2_datamatrices_according_1.R train_mat.csv test_mat.csv  
args <- commandArgs(TRUE);
library("preprocessCore")
library("caret")
#library(affyPLM)

#ref<-read.table("train",header=T,sep=",",row.names=1)
ref<-read.table(args[1],header=T,sep=",",row.names=1)
#test<-read.table("test",header=T,sep=",",row.names=1)
test<-read.table(args[2],header=T,sep=",",row.names=1)
ref1<-as.matrix(ref)
test1<-as.matrix(test)
norm_ref<-normalize.quantiles(ref1)
 colnames(norm_ref) <- colnames(ref1)
 rownames(norm_ref) <- rownames(ref1)
target <- normalize.quantiles.determine.target(norm_ref)
#target
tt <- normalize.quantiles.use.target(test1,target)
#tt
rownames(tt) <- rownames(test1)
colnames(tt) <- colnames(test1)
#tt
#write.table(norm_ref,file="Quant_normalizaed_train.csv",sep=",",quote=F)
write.table(norm_ref,file=args[3],sep=",",quote=F)
#write.table(tt,file="Quant_normalizaed_test.csv",sep=",",quote=F)
write.table(tt,file=args[4],sep=",",quote=F)
#boxplot(as.data.frame(ref1))
#boxplot(as.data.frame(norm_ref))
#boxplot(as.data.frame(test1))
#boxplot(as.data.frame(tt))
 
