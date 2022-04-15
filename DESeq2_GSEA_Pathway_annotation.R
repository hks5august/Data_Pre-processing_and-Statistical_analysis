if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pathview")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")



if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("gage")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("gageData")


library(dplyr)
library(DESeq2)
library("AnnotationDbi")
library("org.Hs.eg.db")
#columns(org.Hs.eg.db)

library(pathview)
library(gage)
library(gageData)


#countDataURL = "https://raw.githubusercontent.com/pine-bio-support/Pathways-GSEA-analysis1/main/GSE37704_sailfish_genecounts.csv"
#countData = read.csv(countDataURL, row.names=1) %>% 
#dplyr::select(-length) %>%  
#as.matrix()

# Import countdata
countData <- read.csv("https://raw.githubusercontent.com/gexijin/iDEP-old-version/master/GSE37704_sailfish_genecounts.csv", row.names=1)
dim(countData)


countData <- as.matrix(countData)
#remove o 0r 1
countData <- countData[rowSums(countData)>1, ]
#dim(countData)
#countData <- as.matrix(countData)

# Import metadata

colData <- read.csv("https://raw.githubusercontent.com/pine-bio-support/Pathways-GSEA-analysis1/main/GSE37704_metadata.csv", row.names=1)
dim(colData)

# Set up the DESeqDataSet Object and run the DESeq pipeline
dds <- DESeqDataSetFromMatrix(countData=countData, colData = colData, design = ~condition)
dds <- DESeq(dds)
dds

# Results
res <- results(dds, contrast=c("condition", "hoxa1_kd", "control_sirna"))
res <- res[order(res$pvalue),]
#res
summary(res)


library("AnnotationDbi")
library("org.Hs.eg.db")
#columns(org.Hs.eg.db)

library(pathview)
library(gage)
library(gageData)

head(res)



res$symbol = mapIds(org.Hs.eg.db,keys=row.names(res), column="SYMBOL",keytype="ENSEMBL",multiVals="first")

res$entrez = mapIds(org.Hs.eg.db,keys=row.names(res), column="ENTREZID",keytype="ENSEMBL",multiVals="first")

res$name =  mapIds(org.Hs.eg.db,keys=row.names(res), column="GENENAME",keytype="ENSEMBL",multiVals="first")



data(kegg.sets.hs)
data(sigmet.idx.hs)

kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]



foldchanges <- res$log2FoldChange
names(foldchanges) <- res_entrez
head(foldchanges)


keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)


attributes(keggres)

lapply(keggres, head)


## Sanity check displaying all pathways data
pathways = data.frame(id=rownames(keggres$greater), keggres$greater)
head(pathways)


pathview(gene.data=foldchanges, pathway.id="hsa04110")


# A different PDF based output of the same data
pathview(gene.data=foldchanges, pathway.id="hsa04110", kegg.native=FALSE)


## Focus on top 5 upregulated pathways here for demo purposes only
keggrespathways <- rownames(keggres$greater)[1:5]

# Extract the IDs part of each string
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids

# Finally, lets pass these IDs in keggresids to the pathview() function to draw plots for all the top 5 pathways.
pathview(gene.data=foldchanges, pathway.id=keggresids, species="hsa")


# GO ontology
data(go.sets.hs)
data(go.subs.hs)
gobpsets = go.sets.hs[go.subs.hs$BP]

gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)

lapply(gobpres, head)

# Reactome Pathway Analysis
#First, Using R, output the list of significant genes at the 0.05 level as a plain text file:
head(res)
 
sig_genes <- res[res$padj <= 0.05 & !is.na(res$padj), "symbol"]

print(paste("Total number of significant genes:", length(sig_genes)))


write.table(sig_genes, file="significant_genes.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
