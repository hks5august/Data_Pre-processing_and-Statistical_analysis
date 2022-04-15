library("Hmisc")
library("lineup")
df<-read.csv("CL184A1_iMEC.csv",header=T, sep=",",row.names=1)

mek<-df[grep("iMEK",rownames(df)),]
egfr<-df[grep("EGFR",rownames(df)),]
full<-df[grep("full",rownames(df)),]
egf<-df[grep("inhEGF",rownames(df)),]
pi3k<-df[grep("PI3K",rownames(df)),]
pkc<-df[grep("PKC",rownames(df)),]

#for(i in seq(from=1, to=37 ,by=1))
#{
#corr<-cor(mek[i],egfr[1:47712,][i])
  
#{write.table(rbind(corr),file="Corr_result_MEK_EGFR.csv",col.names=F,sep = '\t',append = T);#output file
#}
  
#}


######### result by corbet #### 

cor_res<-corbetw2mat(mek, egfr[1:nrow(mek),], what = c("paired"), corthresh = 0.2)
write.table(cor_res,file="Corr_result_MEK_EGFR.csv",col.names=F,sep = '\t',append = T)


