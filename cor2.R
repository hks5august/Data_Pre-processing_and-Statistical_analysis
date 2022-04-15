library(Hmisc)
library("caret")
options(scipen=999)
library(MASS);
args <- commandArgs(TRUE)
pos=read.csv(args[1],check.names=FALSE)
cor1=cor(pos)
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
}
res<-rcorr(as.matrix(pos))
out=flattenCorrMatrix(res$r, res$P)
write.csv(out,file=args[2])
