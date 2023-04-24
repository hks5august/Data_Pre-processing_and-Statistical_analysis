#Venn diagram
#install.packages("DescTools")
#install.packages("UpSetR")
#install.packages("venn")         # Install & load venn package
#install.packages("ggpolypath")   # Install ggpolypath package

library(DescTools)
library(UpSetR)
library("venn")
library("ggpolypath")  
set.seed(7)

setwd("/Users/kaurh8/Documents/CCGA_datasets/CGGA_mRNA_693_samples")

#list_groups()


grade4 <- read.table ("Sig_359_genes_with_sig_padj_logFC_1_5_Grade4.txt", header=T, sep="\t")
grade3 <- read.table ("Sig_2692_genes_with_sig_padj_logFC_1_5_Grade3.txt", header=T, sep="\t")
grade2<- read.table ("Sig_588_genes_with_sig_padj_logFC_1_5_Grade2.txt", header=T, sep="\t")


list_groups<- list(grade4$ID , grade3$ID , grade2$ID )

#?PlotVenn
venn_diag <- PlotVenn(x=list_groups, labels=c("Grade4", "Grade3", "Grade2"), col=SetAlpha(c("blue","red","green"), 0.2))
venn_diag


jpeg(file="Venn_diagram_for_all_grades.jpeg", units="in", width=10, height=10, res=350)
venn_diag
dev.off()


getwd()
