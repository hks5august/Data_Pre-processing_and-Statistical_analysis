options(scipen=999)
library(MASS);
args <- commandArgs(TRUE)
mat<- read.table(args[1], sep=",",header=T)
df<-as.data.frame(mat)
dft <- df[ , grepl( "A" , names( df ) ) ]
dfn <- df[ , grepl( "D" , names( df ) ) ]
write.table(dft,file="tumor",sep="\t");
write.table(dfn,file="normal",sep="\t");
