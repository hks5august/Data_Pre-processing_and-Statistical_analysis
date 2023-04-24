#Venn diagram
install.packages("DescTools")
install.packages("UpSetR")
install.packages("venn")         # Install & load venn package
install.packages("ggpolypath")   # Install ggpolypath package

library(DescTools)
library(UpSetR)
library("venn")
library("ggpolypath")  

setwd("/Users/kaurh8/Documents/CCGA_datasets/mRNA_693_samples/LGG/primary_vs_Recurrent/DGE/Venn_diagram")

List1 <- read.table ("131_genes_list", header=T, sep="\t")
List2 <- read.table ("747_Significant_surv_features_list", header=T, sep="\t")
head(List1)

list_groups<- list(List1$ID, List2$ID)


?PlotVenn
PlotVenn(x=list_groups, labels=c("131 DEGs ", "747 Prognostic"), col=SetAlpha(c("blue","red"), 0.2))


jpeg(file="Venn_16.jpeg", units="in", width=10, height=10, res=350)
PlotVenn(x=list_groups, labels=c("129 DEGs ", "747 Prognostic"), col=SetAlpha(c("blue","red"), 0.2))
dev.off()

