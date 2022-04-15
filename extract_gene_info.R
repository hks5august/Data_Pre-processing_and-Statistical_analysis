library(org.Hs.eg.db)
#library(illuminaHumanv4.db)
#library(hgu133plus2.db)


set.seed(7)
ids <- read.csv(file='sym_ent',header =TRUE, sep = "\t", dec = ".", );

id <-as.character(ids$Entrez_Gene_ID)


#cols1 <- c("ENTREZID", "SYMBOL" , "GENENAME" , "ENSEMBL" , "ENSEMBLPROT" , "ENSEMBLTRANS","UNIPROT", "OMIM" , "REFSEQ" ,"UNIGENE", "UCSCKG" ,  "GO",  "ONTOLOGY", "UCSCKG", "PFAM", "IPI" , "PROSITE","GOALL", "ONTOLOGYALL")


cols1 <- c("ENTREZID", "SYMBOL" , "GENENAME" , "ENSEMBL" , "ENSEMBLPROT" , "ENSEMBLTRANS","UNIPROT", "OMIM" , "REFSEQ" ,"UNIGENE")
cols2 <- c("ENTREZID", "UCSCKG" ,  "GO",  "ONTOLOGY", "UCSCKG", "PFAM", "IPI" , "PROSITE","GOALL", "ONTOLOGYALL")


Ids_gene1 <-select(org.Hs.eg.db, keys=id, columns=cols1, keytype="ENTREZID")

Ids_gene2 <-select(org.Hs.eg.db, keys=id, columns=cols2, keytype="ENTREZID")

final_Ids_gene <- cbind(Ids_gene1,Ids_gene2)

# idsIds_gene<-select(hgu133plus2.db, keys=ENTREZID, columns=c("SYMBOL","GENENAME", "ENTREZID"), keytype="ENTREZID")

#Ids_gene <- select(illuminaHumanv4.db, rownames(data), c("SYMBOL","ENTREZID", "ENSEMBLTRANS", "GENENAME"),"ENTREZID")

write.table(final_Ids_gene, file = "Final_Ids_genes.tsv", quote = FALSE, sep = "\t")
