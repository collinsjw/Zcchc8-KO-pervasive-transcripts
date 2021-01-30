if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
library(edgeR)
library(biomaRt)

#############################
#NIH Zcchc8 KO in SIMS Cells#
#############################

#Get filtered raw read count file 
rc <- read.delim(file = "Filtered/raw/count/file", row.names = 1)

#edgeR on features from SIMS Zcchc8 KO cells
gp <- factor(c(rep("BtKO", 4), rep("WT", 3), rep("ZcKO", 12)))
x <- DGEList(counts = rc, group = gp)
x <- calcNormFactors(x, method = "TMM")
design <- model.matrix(~gp)
x <- estimateCommonDisp(x)
x <- estimateTagwiseDisp(x)
fit <- glmFit(x, design)

#Set n equal to number of features ids; this is to help generate the full list below
n = dim(x$counts)[1]


#Differential features for first group comparison
WvZ_deg <- exactTest(x, pair = c("WT", "ZcKO"))
WvZ <- topTags(WvZ_deg, n = n)
WvZ_res <- as.data.frame(WvZ)


#Differential features for second group comparison
WvB_deg <- exactTest(x, pair = c("WT", "BtKO"))
WvB <- topTags(WvB_deg, n = n)
WvB_res <- as.data.frame(WvB)


#Differential features for third group comparison
BvZ_deg <- exactTest(x, pair = c("BtKO", "ZcKO"))
BvZ <- topTags(BvZ_deg, n = n)
BvZ_res <- as.data.frame(BvZ)

#MDS plots
points <- c(20)
colors <- c("blue", "darkgreen", "red")
color=as.numeric(x$samples$group)
plotMDS(x, col=color, pch=points, cex = 2, pin = c(4, 4), cex.lab=1.5, cex.axis=1.5)


#Annotate ensembl ids with gene names
mus = useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
gene_ids <- getBM(attributes = c("external_gene_name", "ensembl_gene_id", "description"), 
                  filters = "ensembl_gene_id",
                  values = rownames(rc),
                  mart = mus)

WvZ_res <- merge(WvZ_res, gene_ids, by.x = 0, by.y = "ensembl_gene_id")
BvZ_res <- merge(BvZ_res, gene_ids, by.x = 0, by.y = "ensembl_gene_id")
WvB_res <- merge(WvB_res, gene_ids, by.x = 0, by.y = "ensembl_gene_id")


#Write files
write.table(WvZ_res, file = "Your/file/name", sep = "\t", col.names = T, row.names = T, quote = F)

write.table(WvB_res, file = "/Your/file/name", sep = "\t", col.names = T, row.names = T, quote = F)

write.table(BvZ_res, file = "/Your/file/name", sep = "\t", col.names = T, row.names = T, quote = F)

