library(data.table)
library(dplyr)
args <- commandArgs(TRUE);
set.seed(7)
#setwd("/Users/kaurh8/Documents/GDC_TCGA_Biolinks/GDC_All_samples_Data/TCGA-BLCA", row.names=1)
data <- read.table("TCGA-BLCA_Matrix_FPKM_log_data_with_unique_sample_ID_without_SD_tpose", header=T, sep="\t", check.names = F)
#data <- read.table(args[1], header=T, sep="\t", check.names = F)
dim(data)
head(data[2:ncol(data)])
data1 <- data[2:ncol(data)]

#calculate standard deviation in row
row_sd1 <- as.data.frame(apply(data1, 1, sd))                     # Using apply() function
colnames(row_sd1) <- c("SD")
dim(row_sd1)
data2 <- cbind(data, row_sd1)
dim(data2)
head(data2)

#extract ID of samples/genes
#samples <- as.data.frame(data2$gene_name)
samples <- as.data.frame(data2[1])
colnames(samples) <- c("ID")
head(samples)
dim(samples)
dim(unique(samples))
#Extract duplicates
duplicated_samples <- subset(samples,duplicated(samples$ID))
dim(duplicated_samples)
dim(samples)
dim(unique(samples))
head(duplicated_samples)
#uniq dulicated samples
duplicated_samples_uniq <- unique(duplicated_samples)
colnames(duplicated_samples) <- c("ID")
dim(duplicated_samples_uniq)


#create dataframe with sample with duplicated ID
dup_data <- as.data.frame(data2[data2[,1] %in% c(duplicated_samples_uniq$ID), ])
dim(dup_data)
head(dup_data,3)

#create dataframe with sample with non-duplicated ID
#non_dup_data <- as.data.frame(data2[!data2$gene_name %in% c(duplicated_samples_uniq$ID), ])
non_dup_data <- as.data.frame(data2[!data2[,1] %in% c(duplicated_samples_uniq$ID), ])
dim(non_dup_data) 
head(non_dup_data,3)
dim(data2)

#Extract row (data) with max SD values 
#Run loop
dup_SD_max1 = data.frame() #create an empty data frame

for(i in seq(from=1, to=length(t(duplicated_samples_uniq)), by=1)) {
  
  b <- duplicated_samples_uniq[i,1]
  #b <- dup1[3,1]
  print(b)
  
  #extract data  with all rows  with specific gene 
  #i_data <- subset(dup_data, grepl(paste0("SNORA71$", sep=""), dup_data$gene_name))
  #i_data <- subset(dup_data, grepl(paste0("SNORA71$", sep=""), dup_data[,1]))
  i_data <- subset(dup_data, grepl(paste0("^",b,"$", sep=""), dup_data[,1]))

  # Extract row with max SD
  i_data_max <- subset(subset(i_data , SD == max(i_data$SD))) # select rows where SD is maximum
  #i_data_max1 <- subset(subset(i_data_max , SD != 0 )) # remove genes SD is zero 
  # keep first row if multiples rows have same SD
  i_data_max1 <- as.data.frame(i_data_max[rowid(rleid(i_data_max$SD)) == 1L, ])
  head(i_data_max1)
  #head(dup_data,2)
  dup_SD_max1  <- rbind(dup_SD_max1 ,i_data_max1 )
  
  #write.table(i_data_max, file="i_data_max.txt",row.names=F,col.names=F,sep = '\t',append = T, quote = F);#output file
  
}


dim(dup_SD_max1)
dim(non_dup_data)

colnames(non_dup_data) <- colnames(non_dup_data)
colnames(dup_SD_max1) <-  colnames(non_dup_data)
#row.names(non_dup_data) <- non_dup_data$ID
#row.names(dup_SD_max1) <- dup_SD_max1$ID


head(dup_SD_max1,2)


length(dup_SD_max1)

if (length(dup_SD_max1) > 1) {
#final_data_with_SD <- rbind(non_dup_data[2:ncol(non_dup_data)], dup_SD_max1[2:ncol(dup_SD_max1)])
final_data_with_SD <- rbind(non_dup_data, dup_SD_max1)
} else {
  final_data_with_SD <- non_dup_data
}

final_data_without_SD  <- final_data_with_SD[ , -which(names(final_data_with_SD) %in% c("SD"))]
dim(final_data_with_SD)
dim(final_data_without_SD)
head(final_data_with_SD,2)
head(final_data_without_SD,2)

# save files 
#write.table(final_data_with_SD,file="final_data_with_SD.txt",sep="\t",quote=F, row.names=F)

#write.table(final_data_without_SD, ,file="final_data_without_SD.txt",sep="\t",quote=F, row.names=F)
#write.table(cbind("ID"=rownames(final_data_with_SD), final_data_with_SD),file=args[2],sep="\t",quote=F, row.names=F)
#write.table(cbind("ID"=rownames(final_data_without_SD), final_data_without_SD),file=args[3],sep="\t",quote=F, row.names=F)
write.table(final_data_with_SD,file=args[2],sep="\t",quote=F, row.names=F)
write.table(final_data_without_SD,file=args[3],sep="\t",quote=F, row.names=F)

