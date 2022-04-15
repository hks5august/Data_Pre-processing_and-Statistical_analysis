options(scipen=999)
library(MASS);
args <- commandArgs(TRUE)
mat<- read.table(args[1], sep=",",header=T, row.names=1, check.names=F)
mean=apply(mat,1,mean);
sd=apply(mat,1,sd);
write.table(mean,args[2],sep="," , quote=F);
write.table(sd,args[3],sep=",", quote=F);
