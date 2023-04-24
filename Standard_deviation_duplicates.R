library(data.table)
args <- commandArgs(TRUE);
set.seed(7)
setwd("/Users/kaurh8/Documents/GDC_TCGA_Biolinks/GDC_All_samples_Data/TCGA-BLCA", row.names=1)
#data <- read.table("example_data", header=T, sep="\t", check.names = F)
data <- read.table(args[1], header=T, sep="\t", check.names = F)
head(data[2:ncol(data)])
data1 <- data[2:ncol(data)]
row_sd1 <- as.data.frame(apply(data1, 1, sd))                     # Using apply() function
colnames(row_sd1) <- c("SD")
head(row_sd1)
dim(row_sd1)
data2 <- cbind(data, row_sd1)
head(data2)

#extract ID of samples/genes
#samples <- as.data.frame(data2$gene_name)
samples <- as.data.frame(data2[1])
head(samples)
colnames(samples) <- c("ID")
#Extract duplicates
duplicated_samples <- subset(samples,duplicated(samples$ID))

dim(duplicated_samples)
head(duplicated_samples)
duplicated_samples_uniq <- unique(duplicated_samples)
colnames(duplicated_samples) <- c("ID")
dim(duplicated_samples_uniq)


head(data2[,1])

#sel_train <- as.data.frame(data1[,colnames(data1) %in% c(row.names(features)), ])
#dup_data <- as.data.frame(data2[data2$gene_name %in% c(duplicated_samples_uniq$ID), ])
dup_data <- as.data.frame(data2[data2[,1] %in% c(duplicated_samples_uniq$ID), ])
dim(dup_data)
head(dup_data,3)
#non_dup_data <- as.data.frame(data2[!data2$gene_name %in% c(duplicated_samples_uniq$ID), ])
non_dup_data <- as.data.frame(data2[!data2[,1] %in% c(duplicated_samples_uniq$ID), ])
dim(non_dup_data) 
head(non_dup_data,3)


#Run loop
dup_SD_max1 = data.frame() #create a data frame
for(i in seq(from=1, to=length(t(duplicated_samples_uniq)), by=1)) {

b <- duplicated_samples_uniq[i,1]
#b <- dup1[3,1]
print(b)

#extract data  with all rows  with specific gene 
 #i_data <- subset(dup_data, grepl(paste0("SNORA71$", sep=""), dup_data$gene_name))
 #i_data <- subset(dup_data, grepl(paste0("SNORA71$", sep=""), dup_data[,1]))
 i_data <- subset(dup_data, grepl(paste0("^",b,"$", sep=""), dup_data[,1]))
 head(i_data )

# Extract row with max SD
i_data_max <- subset(subset(i_data , SD == max(i_data$SD))) # select rows where SD is maximum
#i_data_max1 <- subset(subset(i_data_max , SD != 0 )) remove genes if 
i_data_max1 <- as.data.frame(i_data_max[rowid(rleid(i_data_max$SD)) == 1L, ])
head(i_data_max1)
#head(dup_data,2)
dup_SD_max1  <- rbind(dup_SD_max1 ,i_data_max1 )

#write.table(i_data_max, file="i_data_max.txt",row.names=F,col.names=F,sep = '\t',append = T, quote = F);#output file

}

dim(dup_SD_max1)
head(dup_SD_max1)

final_data_with_SD <- rbind(non_dup_data, dup_SD_max1)
final_data_without_SD  <- final_data_with_SD[ , -which(names(final_data_with_SD) %in% c("SD"))]
dim(final_data_with_SD)
dim(final_data_without_SD)
head(final_data_with_SD,2)
head(final_data_without_SD,2)

# save files 
#write.table(cbind("ID"=rownames(final_data_with_SD), final_data_with_SD),file="final_data_with_SD.txt",sep="\t",quote=F, row.names=F)
#write.table(cbind("ID"=rownames(final_data_without_SD), final_data_without_SD),file="final_data_without_SD.txt",sep="\t",quote=F, row.names=F)
write.table(cbind("ID"=rownames(final_data_with_SD), final_data_with_SD),file=args[2],sep="\t",quote=F, row.names=F)
write.table(cbind("ID"=rownames(final_data_without_SD), final_data_without_SD),file=file=args[3],sep="\t",quote=F, row.names=F)

