###### remove those columns having missing values ##########
####### Rscript missing.R input output ####
options(scipen=999)
library(MASS);
args <- commandArgs(TRUE)
library("caret")
mat=read.table(args[1],sep=",",check.names=FALSE, header=T, row.names=1)
df<-data.frame(mat)
filteredDescr<-df[sapply(df, function(df) !any(is.na(df)))] 
write.table(filteredDescr,file= args[2],sep=",")
