library(MASS);
args <- commandArgs(TRUE);
mat<- read.csv(args[1], sep=",",header=T,row.names=1);
#matt=(mat+1);
#log_mat=log2(matt)
log_mat=log2(mat)
log_mat<-data.frame(log_mat)
#write.table(log_mat,args[2],sep=",",quote=F);
write.table(data.frame('Gene_ID'=rownames(log_mat), log_mat),file=args[2],sep=',',quote = F,row.names = FALSE)
