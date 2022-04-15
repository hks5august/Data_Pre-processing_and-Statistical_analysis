library(limma);
library(GSVA);
#require("limma");


#sample <-read.table("OUT",sep=",",header=T);
sample <-read.table(args[1],sep=",",header=T);

#exp <-read.table("OUT1",sep=",",header=T,row.names = 1);
exp <-read.table(args[2],sep=",",header=T,row.names = 1);
exp<-as.matrix(exp);

#genes<-read.table("GO_BP1_tpose.csv",header=T,sep="\t",fill=T);
genes<-read.table(args[3],header=T,sep="\t",fill=T);
genesets<-as.list(genes);

##### genearte design model matrix
design <- model.matrix(~ factor(sample$Sample_type))
#design
colnames(design) <- c("ALL", "Non-tumor_vs_HCC")


#gsva(expr, gset.idx.list, annotation, method=c("gsva", "ssgsea", "zscore", "plage"), kcdf=c("Gaussian", "Poisson", "none"), abs.ranking=FALSE, min.sz=1, max.sz=Inf, parallel.sz=0, parallel.type="SOCK", mx.diff=TRUE, tau=switch(method, gsva=1, ssgsea=0.25, NA), ssgsea.norm=TRUE, verbose=TRUE)
gsva_es <- gsva(exp, genesets, method="ssgsea", min.sz=1, max.sz=Inf,mx.diff=1,ssgsea.norm=F);
gsva_es<-data.frame(gsva_es)
write.table(data.frame('ID'=rownames(gsva_es), gsva_es),file=args[4],sep=',',quote = F,row.names = F)
######### Define cut off values for significant enrichment score calculation ####
#adjPvalueCutoff <- 0.01
#logFCcutoff <- log2(2)
adjPvalueCutoff <- args[5];
logFCcutoff <- log2(2);

fit <- lmFit(gsva_es, design)
fit <- eBayes(fit)
#allGeneSets <- topTable(fit, coef="Non-tumor_vs_HCC", number=Inf)
#DEgeneSets <- topTable(fit, coef="Non-tumor_vs_HCC", number=Inf, p.value=adjPvalueCutoff, adjust="BH")
#dim(DEgeneSets) 
#res <- decideTests(fit, p.value=adjPvalueCutoff)
#res
#summary(res)

DEgenes <- topTable(fit, coef="Non-tumor_vs_HCC", number=Inf, p.value=adjPvalueCutoff, adjust="BH", lfc=logFCcutoff)
res <- decideTests(fit, p.value=adjPvalueCutoff, lfc=logFCcutoff)
summary(res)
DEgenes<-data.frame(DEgenes)
write.table(data.frame('ID'=rownames(DEgenes), DEgenes),file=args[6],sep=',',quote = F,row.names = F)
