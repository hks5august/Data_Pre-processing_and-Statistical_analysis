args <- commandArgs(TRUE)
path<- "/Users/harpreetkaur/Documents/Harpreet_projects/GEO/"
## instll libaries ###
#source("https://bioconductor.org/biocLite.R")
#biocLite("GEOquery")
#biocLite("affy")
#biocLite("hgu133plus2.db")
#biocLite("hgu133plus2cdf")
#biocLite("illuminaHumanv4.db")
#biocLite("gcrma")
#biocLite("limma")

####Load libraries
library(GEOquery)
library(affy)
#library(illuminaHumanv4.db)
#library(hgu133plus2.db)
#library(hgu133plus2cdf)
library(oligo)
library(limma)
###### run for multiple files ####

#f<-read.table("files", header=F)
f<-read.table(args[1], header=F)
ff<-as.list(f)
for (i in 1:args[2]) {x <- sapply(ff, "[", i)
# Set working directory for download
setwd(paste(path))
##### download raw data
getGEOSuppFiles(x)
# Set working directory for download
#setwd("/Users/harpreetkaur/Documents/Harpreet_projects/GEO/GSE*")
setwd(paste(path,file = paste0(x),sep=""))
# unpack raw files
zipF <- list.files(pattern="RAW.tar", full.names = T)
untar(zipF, exdir = "data")
## unzip or gunzip files
files = list.files("data/", pattern = NULL)
sapply(paste("data", files, sep = "/"), gunzip)
### set path for data directory
#path2=(paste(path,file = paste0(x),sep=""))
#setwd(paste(path2,file=data,sep="/"))
setwd("data")
wd<-getwd()
cat("Current working dir: ", wd)
files1 <- list.files(pattern = "CEL")
if (length(files1)>1){

cat ("\n CEL file exist \n ")
#if (isTRUE(file.exists(files1))) {

#cat("##### set new directory: ########\n ")
#setwd("data")
#load data
cat("load data")
raw.data = ReadAffy()
celFiles <- list.celfiles()
affyRaw <- read.celfiles(celFiles)

# You might need to install and load a package for the specific array you are using (this example is mouse gene 2.0 ST)
# It may try to load it automatically, but may fail.  Install & load the library manually if this happens.
eset <- rma(affyRaw)

# Finally, save the data to an output file to be used by other programs, etc (Data will be log2 transformed and normalized)
write.exprs(eset,file="oligo_rma_data.txt")

#cat("load data")
#raw.data = ReadAffy()
# perform RMA normalization (log2)
#data_rma.norm = rma(raw.data)
## perform gcrma
#data.gcrma = justGCRMA()
#### or
# data.rma <- justRMA()
##### aletrnatively u can use mass package 	
## data.mass <- mas5(Data)
#raw_data_matrix <- computeExprSet(raw.data, pmcorrect.method="pmonly",summary.method="avgdiff")
############ write data into matrix file
#write.exprs(raw_data_matrix, file="raw_data_matrix.txt")
#write.exprs(data_rma.norm, file="mydata_rma.txt")
#write.exprs(data.gcrma, file="mydata_gcrma.txt")
#write.exprs(data.rma, file="mydata_rma.txt")

#### command to know annotation database name
annot<-annotation(raw.data)
db <-paste0(annot, ".db" )
write.table(db, file = "annoatation_db", quote = FALSE, sep = "\t")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(db,character.only=TRUE)
#source("https://bioconductor.org/biocLite.R")
#biocLite(db,character.only=TRUE)
library(db,character.only=TRUE)
###Extract gene names corresponds to probe IDS
#library(hgu133plus2.db)

#ids_gene1<-select(get(db), keys=probeID, columns=c("SYMBOL", "ENTREZID", "GENENAME"), keytype=regex("PROBEID", ignore_case = TRUE))
Ids_gene1<-select(get(db), rownames(raw.data), columns=c("SYMBOL", "ENTREZID", "GENENAME"), keytype="PROBEID")
#ids_gene1<-select(get(db), keys=PROBEID, columns=c("SYMBOL", "ENTREZID", "TXID", "GENENAME"), keytype="PROBEID")
## Ids_gene<-select(hgu133plus2.db, keys=probeID, columns=c("SYMBOL","GENENAME", "ENTREZID"), keytype="PROBEID")

write.table(Ids_gene1, file = "Ids_gene_affy", quote = FALSE, sep = "\t")



###### for Illumina data ###
#Ids_gene <- select(get(db), rownames(raw.data), c("SYMBOL","ENTREZID", "ENSEMBLTRANS", "GENENAME"),"PROBEID")
#Ids_gene <- select(illuminaHumanv4.db, rownames(data), c("SYMBOL","ENTREZID", "ENSEMBLTRANS", "GENENAME"),"PROBEID")

#write.table(Ids_gene, file = "Ids_gene_illumina", quote = FALSE, sep = "\t")
#or

#Ids_gene <- mapIds(illuminaHumanv4.db, rownames(data), "SYMBOL","PROBEID", multiVals = "list")

}

else {
    cat(".CEL does not exist \n" )
    # Handle this error as appropriate
}
}
