options(scipen=999)
library(MASS);
args <- commandArgs(TRUE)
mat<- read.table(args[1], sep=",",header=T)
matt=t(mat)
write.table(matt,args[2],sep=",");
