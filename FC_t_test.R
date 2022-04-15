set.seed(15)
#options(scipen=999)
args <- commandArgs(TRUE)
#library(matrixStats)

mysample <- read.table(args[1], sep=",", header=TRUE)
gene <- mysample[grep('^Symbol',names(mysample))]
#print(gene);
#t.result <- apply(df[,2:244], 1, function (x) t.test(x[1:152],x[153:243],paired=F))
#p_value <- unlist(lapply(t.result, function(x) x$p.value))



for(i in 1:nrow(mysample))
{

group1 <- mysample[i,grep('^HCC',names(mysample))]
group2 <- mysample[i,grep('^adjacent',names(mysample))]
#gene <- mysample[grep('^Symbol',names(mysample))]
#print(gene);
group1_avg <- sum(group1)/length(group1)
group2_avg <- sum(group2)/length(group2)

gene[i, "Cancer_mean"] <- group1_avg
gene[i, "adjacent_mean"] <- group2_avg
mean_diff <- (group1_avg - group2_avg)

logfc <- abs(log2(group1_avg/group2_avg))


gene[i, "LogFC"] <-  logfc
gene[i, "mean-diff"] <-  mean_diff

temp_data1 <- data.frame(values=c(group1,group2),vars = rep(c("Cancer","adjacent"), times = c(length(group1),length(group2))))
vars1 = c(rep(c("Cancer","adjacent"), times = c(length(group1),length(group2))))
#print (vars1);
values1=c(t(group1),t(group2))
#print (values1);

x <- t.test(values1 ~ vars1, data = temp_data1)
#print(x);
pv <- c(x$p.value)
#print(pv);
fdr <- p.adjust(pv, method = "fdr")
bone <- p.adjust(pv, method = "bonferroni")

#print(fdr);
#print(bone);



y <- wilcox.test(values1 ~ vars1, data = temp_data1)
#print(y);

wilcox <- c(y$p.value)
wilcox_fdr <- p.adjust(wilcox, method = "fdr")
wilcox_bone <- p.adjust(wilcox, method = "bonferroni")


gene[i, "T_test_pval"] <-  pv
gene[i, "T_test_FDR"] <-  fdr
gene[i, "T_test_Bonferroni"] <-  bone

gene[i, "WilcoxTest_pval"] <-  wilcox
gene[i, "wilcox_FDR"] <-  wilcox_fdr
gene[i, "wilcox_Bonferroni"] <-  wilcox_bone


}

write.table(gene,args[2], sep="\t", quote=FALSE, row.name=FALSE)
