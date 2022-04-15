
df<-read.table("tpose",header=T,sep=",")
lb=df[1]
gene=df[7]
ff<-cbind(lb,gene)
normal<-ff[grep("normal", df$Label),]
hcc<-ff[grep("HCC", df$Label),]
hc_freq<-split(hcc, cut(hcc$MARCO, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), include.lowest=TRUE))
normal_freq<-split(normal, cut(normal$MARCO, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), include.lowest=TRUE))
a<-do.call(c, list(hc_freq, normal_freq))
capture.output(a, file = "MARCO.txt")


gene=df[6]
ff<-cbind(lb,gene)
normal<-ff[grep("normal", df$Label),]
hcc<-ff[grep("HCC", df$Label),]
hc_freq<-split(hcc, cut(hcc$MCM3, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), include.lowest=TRUE))
normal_freq<-split(normal, cut(normal$MCM3, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), include.lowest=TRUE))
a<-do.call(c, list(hc_freq, normal_freq))
capture.output(a, file = "MCM3.txt")

gene=df[5]
ff<-cbind(lb,gene)
normal<-ff[grep("normal", df$Label),]
hcc<-ff[grep("HCC", df$Label),]
hc_freq<-split(hcc, cut(hcc$MCM7, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), include.lowest=TRUE))
normal_freq<-split(normal, cut(normal$MCM7, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), include.lowest=TRUE))
a<-do.call(c, list(hc_freq, normal_freq))
capture.output(a, file = "MCM7.txt")


gene=df[10]
ff<-cbind(lb,gene)
normal<-ff[grep("normal", df$Label),]
hcc<-ff[grep("HCC", df$Label),]
hc_freq<-split(hcc, cut(hcc$STEAP3, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), include.lowest=TRUE))
normal_freq<-split(normal, cut(normal$STEAP3, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), include.lowest=TRUE))
a<-do.call(c, list(hc_freq, normal_freq))
capture.output(a, file = "STEAP3.txt")

gene=df[11]
ff<-cbind(lb,gene)
normal<-ff[grep("normal", df$Label),]
hcc<-ff[grep("HCC", df$Label),]
hc_freq<-split(hcc, cut(hcc$ITGA6, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), include.lowest=TRUE))
normal_freq<-split(normal, cut(normal$ITGA6, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), include.lowest=TRUE))
a<-do.call(c, list(hc_freq, normal_freq))
capture.output(a, file = "ITGA6.txt")


gene=df[12]
ff<-cbind(lb,gene)
normal<-ff[grep("normal", df$Label),]
hcc<-ff[grep("HCC", df$Label),]
hc_freq<-split(hcc, cut(hcc$HGFAC, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), include.lowest=TRUE))
normal_freq<-split(normal, cut(normal$HGFAC, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), include.lowest=TRUE))
a<-do.call(c, list(hc_freq, normal_freq))
capture.output(a, file = "HGFAC.txt")


gene=df[14]
ff<-cbind(lb,gene)
normal<-ff[grep("normal", df$Label),]
hcc<-ff[grep("HCC", df$Label),]
hc_freq<-split(hcc, cut(hcc$SSR2, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), include.lowest=TRUE))
normal_freq<-split(normal, cut(normal$SSR2, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), include.lowest=TRUE))
a<-do.call(c, list(hc_freq, normal_freq))
capture.output(a, file = "SSR2.txt")


gene=df[15]
ff<-cbind(lb,gene)
normal<-ff[grep("normal", df$Label),]
hcc<-ff[grep("HCC", df$Label),]
hc_freq<-split(hcc, cut(hcc$STMN1, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), include.lowest=TRUE))
normal_freq<-split(normal, cut(normal$STMN1, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), include.lowest=TRUE))
a<-do.call(c, list(hc_freq, normal_freq))
capture.output(a, file = "STMN1.txt")


gene=df[17]
ff<-cbind(lb,gene)
normal<-ff[grep("normal", df$Label),]
hcc<-ff[grep("HCC", df$Label),]
hc_freq<-split(hcc, cut(hcc$POLD1, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), include.lowest=TRUE))
normal_freq<-split(normal, cut(normal$POLD1, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), include.lowest=TRUE))
a<-do.call(c, list(hc_freq, normal_freq))
capture.output(a, file = "POLD1.txt")


gene=df[21]
ff<-cbind(lb,gene)
normal<-ff[grep("normal", df$Label),]
hcc<-ff[grep("HCC", df$Label),]
hc_freq<-split(hcc, cut(hcc$CXCL12, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), include.lowest=TRUE))
normal_freq<-split(normal, cut(normal$CXCL12, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), include.lowest=TRUE))
a<-do.call(c, list(hc_freq, normal_freq))
capture.output(a, file = "CXCL12.txt")


gene=df[22]
ff<-cbind(lb,gene)
normal<-ff[grep("normal", df$Label),]
hcc<-ff[grep("HCC", df$Label),]
hc_freq<-split(hcc, cut(hcc$ZWINT, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), include.lowest=TRUE))
normal_freq<-split(normal, cut(normal$ZWINT, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), include.lowest=TRUE))
a<-do.call(c, list(hc_freq, normal_freq))
capture.output(a, file = "ZWINT.txt")


gene=df[23]
ff<-cbind(lb,gene)
normal<-ff[grep("normal", df$Label),]
hcc<-ff[grep("HCC", df$Label),]
hc_freq<-split(hcc, cut(hcc$SPATS2, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), include.lowest=TRUE))
normal_freq<-split(normal, cut(normal$SPATS2, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), include.lowest=TRUE))
a<-do.call(c, list(hc_freq, normal_freq))
capture.output(a, file = "SPATS2.txt")


gene=df[24]
ff<-cbind(lb,gene)
normal<-ff[grep("normal", df$Label),]
hcc<-ff[grep("HCC", df$Label),]
hc_freq<-split(hcc, cut(hcc$NSUN5, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), include.lowest=TRUE))
normal_freq<-split(normal, cut(normal$NSUN5, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), include.lowest=TRUE))
a<-do.call(c, list(hc_freq, normal_freq))
capture.output(a, file = "NSUN5.txt")


gene=df[25]
ff<-cbind(lb,gene)
normal<-ff[grep("normal", df$Label),]
hcc<-ff[grep("HCC", df$Label),]
hc_freq<-split(hcc, cut(hcc$MT1E, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), include.lowest=TRUE))
normal_freq<-split(normal, cut(normal$MT1E, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), include.lowest=TRUE))
a<-do.call(c, list(hc_freq, normal_freq))
capture.output(a, file = "MT1E.txt")


gene=df[26]
ff<-cbind(lb,gene)
normal<-ff[grep("normal", df$Label),]
hcc<-ff[grep("HCC", df$Label),]
hc_freq<-split(hcc, cut(hcc$GPSM2, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), include.lowest=TRUE))
normal_freq<-split(normal, cut(normal$GPSM2, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), include.lowest=TRUE))
a<-do.call(c, list(hc_freq, normal_freq))
capture.output(a, file = "GPSM2.txt")


gene=df[27]
ff<-cbind(lb,gene)
normal<-ff[grep("normal", df$Label),]
hcc<-ff[grep("HCC", df$Label),]
hc_freq<-split(hcc, cut(hcc$COL15A1, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), include.lowest=TRUE))
normal_freq<-split(normal, cut(normal$COL15A1, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), include.lowest=TRUE))
a<-do.call(c, list(hc_freq, normal_freq))
capture.output(a, file = "COL15A1.txt")





