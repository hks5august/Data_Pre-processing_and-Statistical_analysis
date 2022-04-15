options(scipen=999)
library(MASS);
args <- commandArgs(TRUE)
mat<- read.table(args[1], sep=",",header=T)
mean=apply(mat,2,mean);
sd=apply(mat,2,sd);
write.matrix(mean,args[2],sep=",");
write.matrix(sd,args[3],sep=",");
