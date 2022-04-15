options(scipen=999);
args <- commandArgs(TRUE);
df <- read.table(args[1], sep=",", header=T, row.names=1, check.names=F);
################ remove rows where all sample have zero value ############
##df2<-df[ !rowSums(df[,colnames(df)[(2:ncol(df))]]==0)==ncol(df)-1, ];
############### remove rows where 0 is present in atleast 220 samples
df2<-df[rowSums(df == 0) <= 220, ]
write.table(df2,args[2],sep="," ,quote=F );
