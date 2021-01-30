library(DESeq2)
library(tidyverse)
library(biomaRt)
library(RColorBrewer)
library(pheatmap)


################################################
# DESeq2 on features from SIMS Zcchc8 KO cells #
################################################

#Get list of raw read files assuming files are in working directory; set working directory as needed
setwd("/Your/raw/read/files/location/")

#Create raw count files if needed.  If you already have a filtered raw count file then start below at line 51.
vm24.file.names <- list.files(getwd(), full.names = F)

#Create list of raw read files to be turned into dataframe
vm24.af <- lapply(vm24.file.names, read.table, sep = "\t", header = F)
vm24.af <- as.data.frame(vm24.af)


#Subset dataframe to create matrix for DESeq2.  My data has 55475 observations and I don't need the first 4 rows.
vm24.af <- vm24.af[5:55475, c(1, seq(4, 76, by = 4))]

#Rename column names
#Start by creating vector of sample names to be used for new column names
wt <- paste0("WT", sprintf("%02.0f", 2:4))
bt <- paste0("Bt", sprintf("%02.0f", 1:4))
zc <- paste0("Zc", sprintf("%02.0f", 1:12))
samples <- c("gene_ID", bt, wt, zc)

#Rename columns.  These are the raw counts.
vm24.af <- rename_at(vm24.af, vars(colnames(vm24.af)), ~samples)
write.table(vm24.af, file = "/Your/file/name", row.names = F, col.names = T, quote = F, sep = "\t")

#Filter dataframe for features with >5 reads in at least one of the samples.
vm24.af <- filter_all(vm24.af, any_vars(. > 5))
write.table(vm24.af, file = "/Your/file/name", row.names = F, col.names = T, quote = F, sep = "\t")


#Get filtered read count file or use the one you just created. 
#Be sure and change the first row of IDs to row the row names or DESeq will not work properly.
rc <- read.delim(file = "/Filtered/raw/count/file", row.names = 1)

#Setup colData for DESeq2
#Prep data frame and groups
groups <- c(rep("BtKO", 4), rep("WT", 3), rep("ZcKO", 12))

setup <- data.frame(ensemble_id = colnames(rc), group = groups, row.names = 1, stringsAsFactors = F)

#Make DESeq2 data set
dds <- DESeqDataSetFromMatrix(countData = rc, colData = setup, design = ~group)

#Run DESeq2
dds <- DESeq(dds)

#Generate normalized counts dataframe
dds <- estimateSizeFactors(dds)
nc <- as.data.frame(counts(dds, normalized = T))

#Fetch results
ZvW_res <- results(dds, contrast = c("group", "ZcKO", "WT"))
ZvW_df <- as.data.frame(ZvW_res, stringsAsFactors = F)

ZvB_res <- results(dds, contrast = c("group", "ZcKO", "BtKO"))
ZvB_df <- as.data.frame(ZvB_res, stringsAsFactors = F)

BvW_res <- results(dds, contrast = c("group", "BtKO", "WT"))
BvW_df <- as.data.frame(BvW_res, stringsAsFactors = F)

#Annotate ensembl ids with gene names
mus = useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
gene_ids <- getBM(attributes = c("external_gene_name", "ensembl_gene_id_version", "description"), 
                     filters = "ensembl_gene_id_version",
                     values = rownames(rc),
                     mart = mus)

ZvW_df <- merge(ZvW_df, gene_ids, by.x = 0, by.y = "ensembl_gene_id_version")
write.table(ZvW_df, file = "/Your/file/name", sep = "\t", col.names = T, row.names = F, quote = F)

ZvB_df <- merge(ZvB_df, gene_ids, by.x = 0, by.y = "ensembl_gene_id_version")
write.table(ZvB_df, file = "/Your/file/name", sep = "\t", col.names = T, row.names = F, quote = F)

BvW_df <- merge(BvW_df, gene_ids, by.x = 0, by.y = "ensembl_gene_id_version")
write.table(BvW_df, file = "/Your/file/name", sep = "\t", col.names = T, row.names = F, quote = F)


######################################## Multi Dimensional Scaling ######################################## 

rld <- varianceStabilizingTransformation(dds, blind=FALSE)
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
mds <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mds, as.data.frame(colData(rld)))
ggplot(mds, aes(X1,X2, color = group)) + geom_point(size=3)



######################################## Heatmaps ######################################## 

#Example of how to make heatmaps.
genes_nc <- as.matrix(nc)

#K-means clustering for heatmaps
clusters <- pheatmap(genes_nc, scale = "row", kmeans_k = 2)
names(clusters$kmeans)
clusterDF <- as.data.frame(factor(clusters$kmeans$cluster))
colnames(clusterDF) <- "Cluster"
OrderByCluster <- genes_nc[order(clusterDF$Cluster), ]

#custom heatmap colors
cust.color <- colorRampPalette(c("navy", "royalblue", "#c5c9c7", "whitesmoke", "firebrick", "red"))(n = 299)

pheatmap(OrderByCluster,
         scale="row", show_rownames = FALSE, cluster_rows = FALSE, color = cust.color)



