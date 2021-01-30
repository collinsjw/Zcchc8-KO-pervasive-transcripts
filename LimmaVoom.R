if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("edgeR", "limma", "statmod"))
library(edgeR)
library(limma)
library(statmod)
library(biomaRt)

#############################
#NIH Zcchc8 KO in SIMS Cells#
#############################

#Get filtered raw read count file.
rc <- read.delim(file = "Filtered/raw/read/count/file", row.names = 1)

gp <- factor(c(rep("BtKO", 4), rep("WT", 3), rep("ZcKO", 12)))
design = model.matrix(~ 0 + gp) 

v1 <- voom(as.matrix(rc), design, normalize="quantile")

#Rename columns in design matrix
colnames(design) <- levels(gp)

#Limma Voom
fit <- lmFit(v1, design)

#Zcchc8 KO vs WT
ZvW_cm <- makeContrasts(ZcKO-WT, levels = design)
ZvW_fit <- contrasts.fit(fit, ZvW_cm)
ZvW_ebays = eBayes(ZvW_fit)

ZvW <- topTable(ZvW_ebays, number = dim(rc)[1])


#Zcchc8 KO vs Btbd7 KO negative control
ZvB_cm <- makeContrasts(ZcKO-BtKO, levels = design)
ZvB_fit <- contrasts.fit(fit, ZvB_cm)
ZvB_ebays = eBayes(ZvB_fit)

ZvB <- topTable(ZvB_ebays, number = dim(rc)[1])


#Btbd7 KO negative control vs WT
BvW_cm <- makeContrasts(BtKO-WT, levels = design)
BvW_fit <- contrasts.fit(fit, BvW_cm)
BvW_ebays = eBayes(BvW_fit)

BvW <- topTable(BvW_ebays, number = dim(rc)[1])



#Annotate ensembl ids with gene names
mus = useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
gene_ids <- getBM(attributes = c("external_gene_name", "ensembl_gene_id", "description"), 
                  filters = "ensembl_gene_id",
                  values = rownames(rc),
                  mart = mus)

ZvW <- merge(ZvW, gene_ids, by.x = 0, by.y = "ensembl_gene_id")
ZvB <- merge(ZvB, gene_ids, by.x = 0, by.y = "ensembl_gene_id")
BvW <- merge(BvW, gene_ids, by.x = 0, by.y = "ensembl_gene_id")


#write files
write.table(ZvW, file = "/Your/file/name", sep = "\t", col.names = T, row.names = T, quote = F)

write.table(ZvB, file = "/Your/file/name", sep = "\t", col.names = T, row.names = T, quote = F)

write.table(BvW, file = "/Your/file/name", sep = "\t", col.names = T, row.names = T, quote = F)

