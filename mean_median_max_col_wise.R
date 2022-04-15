options(scipen=999)
library(MASS);
args <- commandArgs(TRUE)
mat<- read.csv(args[1], sep=",",header=T,row.names=1,check.names=FALSE)
mean=apply(mat,2,mean);
#median=apply(mat,2,median);
#max=apply(mat,2,max);
write.matrix(mean,args[2],sep=",");
#write.matrix(median,args[3],sep=",");
#write.matrix(max,args[4],sep=",");
