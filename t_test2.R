library(matrixStats)
df <- read.csv("final_protein_coding_log_374_can_50normal_gene_symbol_tpose", header=FALSE, sep=",")
t.result <- apply(df[,2:425], 1, function (x) t.test(x[1:374],x[375:424],paired=F))  
df$p_value <- unlist(lapply(t.result, function(x) x$p.value))
df$fdr <- p.adjust(df$p_value, method = "fdr")
df$boneferroni <- p.adjust(df$p_value, method = "bonferroni")

write(df$p_value, file = "pval1",ncolumns = 1, sep = ",")
write(df$fdr, file = "fdr1",ncolumns = 1, sep = ",")
write(df$boneferroni, file = "boneferroni_p_val1",ncolumns = 1, sep = ",")
cancer_mean1 <- rowMeans(df[,2:375])

cancer_sd <- rowSds(as.matrix(df[,2:375]), na.rm=TRUE)

normal_mean1 <- rowMeans(df[,376:425])

normal_sd <- rowSds(as.matrix(df[,376:425]), na.rm=TRUE)

write(cancer_mean1, file = "cancer_mean1",ncolumns = 1, sep = ",")

write(normal_mean1, file = "normal_mean1",ncolumns = 1, sep = ",")

mean_diff_c_n <- (cancer_mean1 - normal_mean1)


write(mean_diff_c_n, file = "mean_diff1",ncolumns = 1, sep = ",")

write(cancer_sd, file = "cancer_sd",ncolumns = 1, sep = ",")

write(normal_sd, file = "normal_sd",ncolumns = 1, sep = ",")


