options(scipen=999)
args <- commandArgs(TRUE)
library("OmicsMarkeR")
dat <- readLines(args[1])
dat <- strsplit(dat, ",")
mat=as.matrix(dat)
max.length <- max(sapply(dat, length))
l <- lapply(dat, function(v) { c(v, rep(NA, max.length-length(v)))})
a=do.call(cbind, l)
param <- paste0("V",1:10)
param1 <- paste0("V",1:10)
for (i in 1:10)
        {
                for (j in 1:10)
                        {aa=jaccard(a[,i], a[,j]);cat(paste(aa,","))
}
cat("\n");
        }
