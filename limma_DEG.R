library("limma")

##### design contains lebels mainly, where labels in the 2nd column and 1st column contain all samples (like 1)#######

design <- model.matrix(~ factor(leukemia_es$subtype))


adjPvalueCutoff <- 0.01

logFCcutoff <- log2(1.5)

colnames(design) <- c("ALL", "MLLvsALL")


######### apply limma for DEG analysis ######
fit <- lmFit(leukemia_filtered_eset, design)

fit <- eBayes(fit)



allGenes <- topTable(fit, coef="MLLvsALL", number=Inf)

DEgenes <- topTable(fit, coef="MLLvsALL", number=Inf, p.value=adjPvalueCutoff, adjust="BH", lfc=logFCcutoff)

res <- decideTests(fit, p.value=adjPvalueCutoff, lfc=logFCcutoff)



