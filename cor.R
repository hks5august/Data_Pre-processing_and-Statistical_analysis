options(scipen=999)
library(MASS);
args <- commandArgs(TRUE)
mat<- read.csv(args[1],header=T)
cormat=cor(mat)
write.csv(cormat,args[2])
