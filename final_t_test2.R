options(scipen=999)
args <- commandArgs(TRUE)
library(matrixStats)
df <- read.csv(args[1], header=T, sep="\t")
t.result <- apply(df[,2:244], 1, function (x) t.test(x[1:152],x[153:243],paired=F))  
df$p_value <- unlist(lapply(t.result, function(x) x$p.value))
df$fdr <- p.adjust(df$p_value, method = "fdr")
df$boneferroni <- p.adjust(df$p_value, method = "bonferroni")

write.table(df$p_value, file = "pval1",col.names = c("p-value"),row.names=F,quote=F)
write.table(df$fdr, file = "fdr1",col.names = c("FDR"),row.names=F,quote=F)
write.table(df$boneferroni, file = "boneferroni_p_val1",col.names = c("Bonferroni_p-value"),row.names=F,quote=F)
cancer_mean1 <- rowMeans(df[,2:153])

cancer_sd <- rowSds(as.matrix(df[,2:153]), na.rm=TRUE)

normal_mean1 <- rowMeans(df[,154:244])

normal_sd <- rowSds(as.matrix(df[,154:144]), na.rm=TRUE)

write.table(cancer_mean1, file = "cancer_mean1",col.names = c("Mean_in_Cancer"),row.names=F,quote=F)

write.table(normal_mean1, file = "normal_mean1",col.names = c("Mean_in_Normal"),row.names=F,quote=F)

mean_diff_c_n <- (cancer_mean1 - normal_mean1)


write.table(mean_diff_c_n, file = "mean_diff1",col.names = c("Mean_diff"),row.names=F,quote=F)

write.table(cancer_sd, file = "cancer_sd",col.names = c("Cancer_SD"),row.names=F,quote=F)

write.table(normal_sd, file = "normal_sd",col.names = c("Normal_SD"),row.names=F,quote=F)
