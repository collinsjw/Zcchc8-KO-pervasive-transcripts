library(tidyverse)
library(venn)

##########################################################################################
##########################Predicted Genes Zcchc8 vs WT####################################
##########################################################################################
ngd <- read.delim("GSE165689/genes/ZcKOvsWT_genes_deseq2.txt")
nge <- read.delim("GSE165689/genes/ZcKOvsWT_genes_edgeR.txt")
ngl <- read.delim("GSE165689/genes/ZcKOvsWT_genes_limmaVoom.txt")

ngd_pos <- filter(ngd, ngd$log2FoldChange > 0.58 & ngd$padj <= 0.01)
ngd_neg <- filter(ngd, ngd$log2FoldChange < -0.58 & ngd$padj <= 0.01)

nge_pos <- filter(nge, nge$logFC > 0.58 & nge$FDR <= 0.01)
nge_neg <- filter(nge, nge$logFC < -0.58 & nge$FDR <= 0.01)

ngl_pos <- filter(ngl, ngl$logFC > 0.58 & ngl$adj.P.Val <= 0.01)
ngl_neg <- filter(ngl, ngl$logFC < -0.58 & ngl$adj.P.Val <= 0.01)

g <- ngd_pos$external_gene_name
h <- nge_pos$external_gene_name
i <- ngl_pos$external_gene_name

j <- ngd_neg$external_gene_name
k <- nge_neg$external_gene_name
l <- ngl_neg$external_gene_name

nih.pos.genes <- list("DESeq2" = g, "edgeR" = h, "limmaVoom" = i)
npg.venn <- venn(nih.pos.genes, zcolor = "style")

npg.int <- attributes(npg.venn)
npg.int <- npg.int[["intersections"]][[6]]

gene_vector <- scan(file = "GSE165689/genes/CommonNihGenes.txt", what = "character")
gene_vector <- paste0(gene_vector, "\\b")
genes <- filter(ngd_pos, grepl(paste(gene_vector, collapse = "|"), ngd_pos$external_gene_name))
pred.genes <- filter(genes, grepl("Gm", genes$external_gene_name) | grepl("Rik", genes$external_gene_name) |
                      grepl("expressed sequence", genes$description) | grepl("cDNA sequence", genes$description) |
                      grepl("DNA segment", genes$description))

nih.neg.genes <- list("DESeq2" = j, "edgeR" = k, "limmaVoom" = l)
nng.venn <- venn(nih.neg.genes, zcolor = "style")
nng.int <- attributes(nng.venn)
nng.int <- as.data.frame(nng.int[["intersections"]][[5]])

##########################################################################################
##########################Predicted Genes Btbd7 vs WT#####################################
##########################################################################################
ngd <- read.delim("GSE165689/genes/BtKOvsWT_genes_deseq2.txt")
nge <- read.delim("GSE165689/genes/BtKOvsWT_genes_edgeR.txt")
ngl <- read.delim("GSE165689/genes/BtKOvWT_genes_limmaVoom.txt")

ngd_pos <- filter(ngd, ngd$log2FoldChange > 0.58 & ngd$padj <= 0.01)
ngd_neg <- filter(ngd, ngd$log2FoldChange < -0.58 & ngd$padj <= 0.01)

nge_pos <- filter(nge, nge$logFC > 0.58 & nge$FDR <= 0.01)
nge_neg <- filter(nge, nge$logFC < -0.58 & nge$FDR <= 0.01)

ngl_pos <- filter(ngl, ngl$logFC > 0.58 & ngl$adj.P.Val <= 0.01)
ngl_neg <- filter(ngl, ngl$logFC < -0.58 & ngl$adj.P.Val <= 0.01)

g <- ngd_pos$external_gene_name
h <- nge_pos$external_gene_name
i <- ngl_pos$external_gene_name

j <- ngd_neg$external_gene_name
k <- nge_neg$external_gene_name
l <- ngl_neg$external_gene_name

nih.pos.genes <- list("DESeq2" = g, "edgeR" = h, "limmaVoom" = i)
npg.venn <- venn(nih.pos.genes, zcolor = "style")

npg.int <- attributes(npg.venn)
npg.int <- as.data.frame(npg.int[["intersections"]][[6]])

nih.neg.genes <- list("DESeq2" = j, "edgeR" = k, "limmaVoom" = l)
nng.venn <- venn(nih.neg.genes, zcolor = "style")
nng.int <- attributes(nng.venn)
nng.int <- as.data.frame(nng.int[["intersections"]][[7]])
