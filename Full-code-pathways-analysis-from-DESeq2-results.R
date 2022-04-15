#Load DESeq2 results file as input 
result <- read.table("https://raw.githubusercontent.com/pine-bio-support/Pathways-GSEA-analysis1/main/Deseq2_results_readcount_new_DeSeq2_SignifGene_by_pvalue.txt", sep="\t", header = TRUE, row.names=1)
#Check top 10 lines in results 
head(result)

#Load libraries
library(gage)
library(gageData)
library(pathview)

#list of KEGG pathways genes
data(kegg.sets.hs)

#Index of numbers of pathway
data(sigmet.idx.hs)

#Extract  cleaner gene set of important pathways
kegg.sets.hs <- kegg.sets.hs[sigmet.idx.hs]

#Extract the fold changes from the results
foldchanges <- result$log2FoldChange

#Check by head
head(foldchanges)

#Extract gene ids
names(foldchanges) <- result$ENTREZID

#Check top 10 values
head(foldchanges)

#Get all KEGG pathways associated with genes
kegg_result <- gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)

#Write results of all pathways in a file
write.table(kegg_result, file="all_pathways.txt", sep="\t", row.names=FALSE,quote=FALSE)

#Check attributes
attributes(kegg_result)

#Check upregulated pathways
head(kegg_result$greater)

#check down regulated pathways
head(kegg_result$less)

#Extract all upregulated pathways data
Upreg_pathways <- data.frame(Pathway_id=rownames(kegg_result$greater), kegg_result$greater)

#write all upregulated pathways into a file
write.table(Upreg_pathways, file="upregulated_pathways.txt", sep="\t", row.names=F, quote=FALSE)

#Extract all downregulated pathways data
Downreg_pathways <- data.frame(Pathway_id=rownames(kegg_result$less), kegg_result$less)

#write all downregulated pathways into a file
write.table(Downreg_pathways, file="Downregulated_pathways.txt", sep="\t", row.names=F, quote=FALSE)

#Visualize one upregulated pathway
pathview(gene.data=foldchanges, pathway.id="hsa04514")

#Obtain a different PDF based output of the same data
pathview(gene.data=foldchanges, pathway.id="hsa04514", kegg.native=FALSE)

#Extract the top 10 upregulated pathways 
keggrespathways_up <- rownames(kegg_result$greater)[1:10]

#Extract the IDs part of each string
keggresids_up <- substr(keggrespathways_up, start=1, stop=8)
head(keggresids_up)

#write Ids in a file
write.table(keggresids_up, file="top10_uppathways_ids.txt", sep="\t", row.names=F, quote=FALSE)

#Draw plots for top 10 upregulated pathways
pathview(gene.data=foldchanges, pathway.id=keggresids_up, species="hsa")

#Extract top 10 down regulated pathways 
keggrespathways_down <- rownames(kegg_result$less)[1:10]

#Extract the IDs 
keggresids_down = substr(keggrespathways_down, start=1, stop=8)
head(keggresids_down)

#write Ids in a file
write.table(keggresids_down, file="top10_uppathways_ids.txt", sep="\t", row.names=F, quote=FALSE)

#Draw plots for top 10 downregulated pathways
pathview(gene.data=foldchanges, pathway.id=keggresids_down, species="hsa")

#GO Biological Processes
data(go.sets.hs)
data(go.subs.hs)
go_BP_sets <- go.sets.hs[go.subs.hs$BP]
go_BP_res <- gage(foldchanges, gsets=go_BP_sets, same.dir=TRUE)
lapply(go_BP_res, head)

#write all GO Biological processes (BP) into a file
write.table(go_BP_res, file="GO_BP_result.txt", sep="\t", row.names=FALSE, quote=FALSE)
Upregulated GO BP

#Extract all upregulated GO Biological processes
Upreg_GO_BP <- data.frame(GO_id=rownames(go_BP_res$greater), go_BP_res$greater)

#write all upregulated GO BP into a file
write.table(Upreg_GO_BP, file="Upreg_GO_BP.txt", sep="\t", row.names=F, quote=FALSE)

#Downregulated GO BP
#Extract all downregulated GO Biological processes
Downreg_GO_BP <- data.frame(GO_id=rownames(go_BP_res$less), go_BP_res$less)

#write all upregulated GO BP into a file
write.table(Downreg_GO_BP, file="Downreg_GO_BP.txt", sep="\t", row.names=F, quote=FALSE)


# GO Molecular functions (MF)
data(go.sets.hs)
data(go.subs.hs)
go_MF_sets = go.sets.hs[go.subs.hs$MF]
go_MF_res = gage(foldchanges, gsets=go_MF_sets, same.dir=TRUE)
lapply(go_MF_res, head)

#write all GO Molecular functions (MF) into a file
write.table(go_MF_res, file="GO_MF_result.txt", sep="\t", row.names=FALSE, quote=FALSE)

#Extract all upregulated GO Molecular functions
Upreg_GO_MF <- data.frame(GO_id=rownames(go_MF_res$greater), go_MF_res$greater)

#write all upregulated GO MF into a file
write.table(Upreg_GO_MF, file="Upreg_GO_MF.txt", sep="\t", row.names=F, quote=FALSE)

#Extract all downregulated GO molecular functions
Downreg_GO_MF <- data.frame(GO_id=rownames(go_MF_res$less), go_MF_res$less)

#write all upregulated GO MF into a file
write.table(Downreg_GO_MF, file="Downreg_GO_MF.txt", sep="\t", row.names=F, quote=FALSE)
Cellular Components (CC)

#GO cellular components (CC)
data(go.sets.hs)
data(go.subs.hs)
go_CC_sets = go.sets.hs[go.subs.hs$CC]
go_CC_res = gage(foldchanges, gsets=go_CC_sets, same.dir=TRUE)
lapply(go_CC_res, head)

#write all GO Cellular components (CC) into a file
write.table(go_CC_res, file="GO_CC_result.txt", sep="\t", row.names=FALSE, quote=FALSE)

#Extract all upregulated GO cellular components
Upreg_GO_CC <- data.frame(GO_id=rownames(go_CC_res$greater), go_CC_res$greater)

#write all upregulated GO CC into a file
write.table(Upreg_GO_CC, file="Upreg_GO_CC.txt", sep="\t", row.names=F, quote=FALSE)

#Extract all downregulated GO cellular components
Downreg_GO_CC <- data.frame(GO_id=rownames(go_CC_res$less), go_CC_res$less)

#Write all upregulated GO CC into a file
write.table(Downreg_GO_CC, file="Downreg_GO_CC.txt", sep="\t", row.names=F, quote=FALSE)




