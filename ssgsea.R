#!/usr/bin/env Rscript
#install.packages("devtools")
#library(devtools)
#install_github("rcastelo/GSVA")
library(limma);
library(GSVA);
genes<-read.table("gene_list",header=T,sep="\t")

exp<-read.table("mat1.csv",sep=",",header=T,row.names = 1)

exp<-as.matrix(exp)
genesets<-as.list(genes)

GO_BP<-read.table("GO_BP.csv",sep=",",header=T,row.names = 1)
GO_MF<-read.table("GO_MF.csv",sep=",",header=T,row.names = 1)
GO_CC<-read.table("GO_CC.csv",sep=",",header=T,row.names = 1)


GO_BP<-as.list(GO_BP)
GO_CC<-as.list(GO_CC)
GO_MF<-as.list(GO_MF)

#gsva(expr, gset.idx.list, annotation, method=c("gsva", "ssgsea", "zscore", "plage"), kcdf=c("Gaussian", "Poisson", "none"), abs.ranking=FALSE, min.sz=1, max.sz=Inf, parallel.sz=0, parallel.type="SOCK", mx.diff=TRUE, tau=switch(method, gsva=1, ssgsea=0.25, NA), ssgsea.norm=TRUE, verbose=TRUE)
gsva_es_genes <- gsva(exp, genesets, method="ssgsea", min.sz=1, max.sz=Inf,mx.diff=1)

gsva_es_GO_BP <- gsva(exp, GO_BP, method="ssgsea", min.sz=1, max.sz=Inf,mx.diff=1)

gsva_es_GO_MF <- gsva(exp, GO_MF, method="ssgsea", min.sz=1, max.sz=Inf,mx.diff=1)

gsva_es_GO_CC <- gsva(exp, GO_CC, method="ssgsea", min.sz=1, max.sz=Inf,mx.diff=1)

#write.table(gsva_es,args[4], sep="\t", quote=FALSE, row.name=FALSE)
write.table(gsva_es,gsva_ES_score.tsv, sep="\t", quote=FALSE, row.name=T)

write.table(gsva_es_GO_BP,gsva_ES_score_GO_BP.tsv, sep="\t", quote=FALSE, row.name=T)

write.table(gsva_es_GO_CC,gsva_ES_score_GO_CC.tsv, sep="\t", quote=FALSE, row.name=T)

write.table(gsva_es_GO_MF,gsva_ES_score_GO_MF.tsv, sep="\t", quote=FALSE, row.name=T)

