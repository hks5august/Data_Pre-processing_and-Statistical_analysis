options(scipen=999)
library(MASS);
args <- commandArgs(TRUE)
a=read.csv(args[1],header=T,sep=",")
agrregate=aggregate(. ~ id, data = a, mean)
agrregate<-data.frame(agrregate)
write.table(data.frame('Gene_ID'=rownames(agrregate), agrregate),file=args[2],sep=',',quote = F,row.names = FALSE)
#write.table(agrregate,file=args[2],sep=",",quote=F)
