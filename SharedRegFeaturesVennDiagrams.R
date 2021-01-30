library(tidyverse)
library(venn)
library(ggpolypath)

setwd("/Your/differential/expression/results")

ncd <- read.delim("GSE165689/ZcKOvsWT_ctcf_deseq2.txt")
nce <- read.delim("GSE165689/ZcKOvsWT_ctcf_edgeR.txt")
ncl <- read.delim("GSE165689/ZcKOvsWT_ctcf_limmaVoom.txt")

hcd <- read.delim("GSE126108/ctcf_deseq2.txt")
hce <- read.delim("GSE126108/ctcf_edgeR.txt")
hcl <- read.delim("GSE126108/ctcf_limmaVoom.txt")

scd <- read.delim("GSE127790/ctcf_deseq2.txt")
sce <- read.delim("GSE127790/ctcf_edgeR.txt")
scl <- read.delim("GSE127790/ctcf_limmaVoom.txt")

ncd <- rownames_to_column(ncd, "ID")
nce <- rownames_to_column(nce, "ID")
ncl <- rownames_to_column(ncl, "ID")

hcd <- rownames_to_column(hcd, "ID")
hce <- rownames_to_column(hce, "ID")
hcl <- rownames_to_column(hcl, "ID")

scd <- rownames_to_column(scd, "ID")
sce <- rownames_to_column(sce, "ID")
scl <- rownames_to_column(scl, "ID")



#################################################
# Filter sets based on log2FC > 0.58 or < -0.58 #
#################################################
#GSE165689 Zcchc8 KO SIMS Cells

ncd_pos <- filter(ncd, ncd$log2FoldChange > 0.58 & ncd$padj <= 0.01)
ncd_neg <- filter(ncd, ncd$log2FoldChange < -0.58 & ncd$padj <= 0.01)

nce_pos <- filter(nce, nce$logFC > 0.58 & nce$FDR <= 0.01)
nce_neg <- filter(nce, nce$logFC < -0.58 & nce$FDR <= 0.01)

ncl_pos <- filter(ncl, ncl$logFC > 0.58 & ncl$adj.P.Val <= 0.01)
ncl_neg <- filter(ncl, ncl$logFC < -0.58 & ncl$adj.P.Val <= 0.01)


#GSE126108 KO mouse
hcd_pos <- filter(hcd, hcd$log2FoldChange > 0.58 & hcd$padj <= 0.01)
hcd_neg <- filter(hcd, hcd$log2FoldChange < -0.58 & hcd$padj <= 0.01)

hce_pos <- filter(hce, hce$logFC > 0.58 & hce$FDR <= 0.01)
hce_neg <- filter(hce, hce$logFC < -0.58 & hce$FDR <= 0.01)

hcl_pos <- filter(hcl, hcl$logFC > 0.58 & hcl$adj.P.Val <= 0.01)
hcl_neg <- filter(hcl, hcl$logFC < -0.58 & hcl$adj.P.Val <= 0.01)


#GSE127790 KO ES cells
scd_pos <- filter(scd, scd$log2FoldChange > 0.58 & scd$padj <= 0.01)
scd_neg <- filter(scd, scd$log2FoldChange < -0.58 & scd$padj <= 0.01)

sce_pos <- filter(sce, sce$logFC > 0.58 & sce$FDR <= 0.01)
sce_neg <- filter(sce, sce$logFC < -0.58 & sce$FDR <= 0.01)

scl_pos <- filter(scl, scl$logFC > 0.58 & scl$adj.P.Val <= 0.01)
scl_neg <- filter(scl, scl$logFC < -0.58 & scl$adj.P.Val <= 0.01)


######################################
# GSE165689 SIMS Cells Venn Diagrams #
######################################

a <- ncd_pos$ID
b <- nce_pos$ID
c <- ncl_pos$ID

d <- ncd_neg$ID
e <- nce_neg$ID
f <- ncl_neg$ID


#GSE165689 Positive ctcf
GSE165689.pos.ctcf <- list("DESeq2" = a, "edgeR" = b, "Limma-Voom" = c)
npc.venn <- venn(GSE165689.pos.ctcf, zcolor = "style", ilcs = 1, sncs = 1.2)

npc.int <- attributes(npc.venn)
npc.int <- npc.int[["intersections"]][[7]]


#GSE165689 Negative ctcf
GSE165689.neg.ctcf <- list("DESeq2" = d, "edgeR" = e, "limmaVoom" = f)
nnc.venn <- venn(GSE165689.neg.ctcf, zcolor = "style")

nnc.int <- attributes(nnc.venn)
nnc.int <- nnc.int[["intersections"]][[6]]


###################################
# GSE126108 KO MouseVenn Diagrams #
###################################

hop.a <- hcd_pos$ID
hop.b <- hce_pos$ID
hop.c <- hcl_pos$ID

hop.d <- hcd_neg$ID
hop.e <- hce_neg$ID
hop.f <- hcl_neg$ID

#GSE126108 Positive ctcf
hop.pos.ctcf <- list("DESeq2" = hop.a, "edgeR" = hop.b, "limmaVoom" = hop.c)
hpc.venn <- venn(hop.pos.ctcf, zcolor = "style")

hpc.int <- attributes(hpc.venn)
hpc.int <- hpc.int[["intersections"]][[7]]

#GSE126108 Negative ctcf
hop.neg.ctcf <- list("DESeq2" = hop.d, "edgeR" = hop.e, "limmaVoom" = hop.f)
hnc.venn <- venn(hop.neg.ctcf, zcolor = "style")

hnc.int <- attributes(hnc.venn)
hnc.int <- hnc.int[["intersections"]][[5]]

######################################
#GSE127790 KO ES Cells Venn Diagrams #
######################################

sh.a <- scd_pos$ID
sh.b <- sce_pos$ID
sh.c <- scl_pos$ID

sh.d <- scd_neg$ID
sh.e <- sce_neg$ID
sh.f <- scl_neg$ID


#GSE127790 Positive ctcf
sh.pos.ctcf <- list("DESeq2" = sh.a, "edgeR" = sh.b, "limmaVoom" = sh.c)
spc.venn <- venn(sh.pos.ctcf, zcolor = "style")

spc.int <- attributes(spc.venn)
spc.int <- spc.int[["intersections"]][[6]]

#GSE127790 Negative ctcf
sh.neg.ctcf <- list("DESeq2" = sh.d, "edgeR" = sh.e, "limmaVoom" = sh.f)
snc.venn <- venn(sh.neg.ctcf, zcolor = "style")

snc.int <- attributes(snc.venn)
snc.int <- snc.int[["intersections"]][[5]]



###################################################
# Overlap between GSE165689 GSE126108 and GSE127790   
###################################################
#CTCF
#########################################################################################################################################################
#Positive
#Venn diagram for GSE165689 and GSE126108 CTCF
nh.ctcf <- list("GSE165689 ctcf" = npc.int, "GSE126108 ctcf" = hpc.int)
nh.ctcf.venn <- venn(nh.ctcf, zcolor = "style")
nh.ctcf.int <- attributes(nh.ctcf.venn)
nh.ctcf.int <- nh.ctcf.int[["intersections"]][[3]]

#Venn diagram for GSE165689 and GSE127790 ctcf
ns.ctcf <- list("GSE165689 ctcf" = npc.int, "GSE127790 ctcf" = spc.int)
ns.ctcf.venn <- venn(ns.ctcf, zcolor = "style")
ns.ctcf.int <- attributes(ns.ctcf.venn)
ns.ctcf.int <- ns.ctcf.int[["intersections"]][[3]]

#Venn diagram for GSE165689, GSE126108, and GSE127790 ctcf
nhs.ctcf <- list("GSE165689 ctcf" = npc.int, "GSE127790 ctcf" = spc.int, "GSE126108 ctcf" = hpc.int)
nhs.ctcf.venn <- venn(nhs.ctcf, zcolor = "style", ilcs = 1)
nhs.ctcf.int <- attributes(nhs.ctcf.venn)
nhs.ctcf.int <- nhs.ctcf.int[["intersections"]][[5]]

#Negative
#Venn diagram for GSE165689 and GSE126108 ctcf
nh.ctcf <- list("GSE165689 ctcf" = nnc.int, "GSE126108 ctcf" = hnc.int)
nh.ctcf.venn <- venn(nh.ctcf, zcolor = "style")
nh.ctcf.int <- attributes(nh.ctcf.venn)
nh.ctcf.int <- nh.ctcf.int[["intersections"]][[3]]

#Venn diagram for GSE165689 and GSE127790 ctcf
ns.ctcf <- list("GSE165689 ctcf" = nnc.int, "GSE127790 ctcf" = snc.int)
ns.ctcf.venn <- venn(ns.ctcf, zcolor = "style")
ns.ctcf.int <- attributes(ns.ctcf.venn)
ns.ctcf.int <- ns.ctcf.int[["intersections"]][[3]]

#Venn diagram for GSE165689, GSE126108, and GSE127790 ctcf
nhs.ctcf <- list("GSE165689 ctcf" = nnc.int, "GSE127790 ctcf" = snc.int, "GSE126108 ctcf" = hnc.int)
nhs.ctcf.venn <- venn(nhs.ctcf, zcolor = "style", ilcs = 1)
nhs.ctcf.int <- attributes(nhs.ctcf.venn)
nhs.ctcf.int <- nhs.ctcf.int[["intersections"]][[]]


#enhancers
#########################################################################################################################################################
ncd <- read.delim("GSE165689/ZcKOvsWT_enhancers_deseq2.txt")
nce <- read.delim("GSE165689/ZcKOvsWT_enhancers_edgeR.txt")
ncl <- read.delim("GSE165689/ZcKOvsWT_enhancers_limmaVoom.txt")

hcd <- read.delim("GSE126108 /enhancers_deseq2.txt")
hce <- read.delim("GSE126108 /enhancers_edgeR.txt")
hcl <- read.delim("GSE126108 /enhancers_limmaVoom.txt")

scd <- read.delim("GSE127790/enhancers_deseq2.txt")
sce <- read.delim("GSE127790/enhancers_edgeR.txt")
scl <- read.delim("GSE127790/enhancers_limmaVoom.txt")

ncd <- rownames_to_column(ncd, "ID")
nce <- rownames_to_column(nce, "ID")
ncl <- rownames_to_column(ncl, "ID")

hcd <- rownames_to_column(hcd, "ID")
hce <- rownames_to_column(hce, "ID")
hcl <- rownames_to_column(hcl, "ID")

scd <- rownames_to_column(scd, "ID")
sce <- rownames_to_column(sce, "ID")
scl <- rownames_to_column(scl, "ID")



#################################################
# Filter sets based on log2FC > 0.58 or < -0.58 #
#################################################
#GSE165689 Zcchc8 KO SIMS Cells

ncd_pos <- filter(ncd, ncd$log2FoldChange > 0.58 & ncd$padj <= 0.01)
ncd_neg <- filter(ncd, ncd$log2FoldChange < -0.58 & ncd$padj <= 0.01)

nce_pos <- filter(nce, nce$logFC > 0.58 & nce$FDR <= 0.01)
nce_neg <- filter(nce, nce$logFC < -0.58 & nce$FDR <= 0.01)

ncl_pos <- filter(ncl, ncl$logFC > 0.58 & ncl$adj.P.Val <= 0.01)
ncl_neg <- filter(ncl, ncl$logFC < -0.58 & ncl$adj.P.Val <= 0.01)


#GSE126108 KO mouse
hcd_pos <- filter(hcd, hcd$log2FoldChange > 0.58 & hcd$padj <= 0.01)
hcd_neg <- filter(hcd, hcd$log2FoldChange < -0.58 & hcd$padj <= 0.01)

hce_pos <- filter(hce, hce$logFC > 0.58 & hce$FDR <= 0.01)
hce_neg <- filter(hce, hce$logFC < -0.58 & hce$FDR <= 0.01)

hcl_pos <- filter(hcl, hcl$logFC > 0.58 & hcl$adj.P.Val <= 0.01)
hcl_neg <- filter(hcl, hcl$logFC < -0.58 & hcl$adj.P.Val <= 0.01)


#GSE127790 KO ES cells
scd_pos <- filter(scd, scd$log2FoldChange > 0.58 & scd$padj <= 0.01)
scd_neg <- filter(scd, scd$log2FoldChange < -0.58 & scd$padj <= 0.01)

sce_pos <- filter(sce, sce$logFC > 0.58 & sce$FDR <= 0.01)
sce_neg <- filter(sce, sce$logFC < -0.58 & sce$FDR <= 0.01)

scl_pos <- filter(scl, scl$logFC > 0.58 & scl$adj.P.Val <= 0.01)
scl_neg <- filter(scl, scl$logFC < -0.58 & scl$adj.P.Val <= 0.01)


################################
# GSE165689 SIMS Cells Venn Diagrams #
################################

a <- ncd_pos$ID
b <- nce_pos$ID
c <- ncl_pos$ID

d <- ncd_neg$ID
e <- nce_neg$ID
f <- ncl_neg$ID


#GSE165689 Positive enhancers
GSE165689.pos.enhancers <- list("DESeq2" = a, "edgeR" = b, "Limma-Voom" = c)
npc.venn <- venn(GSE165689.pos.enhancers, zcolor = "style", ilcs = 1, sncs = 1.2)

npc.int <- attributes(npc.venn)
npc.int <- npc.int[["intersections"]][[7]]


#GSE165689 Negative enhancers
GSE165689.neg.enhancers <- list("DESeq2" = d, "edgeR" = e, "limmaVoom" = f)
nnc.venn <- venn(GSE165689.neg.enhancers, zcolor = "style")

nnc.int <- attributes(nnc.venn)
nnc.int <- nnc.int[["intersections"]][[5]]


#################################
# GSE126108 KO MouseVenn Diagrams #
#################################

hop.a <- hcd_pos$ID
hop.b <- hce_pos$ID
hop.c <- hcl_pos$ID

hop.d <- hcd_neg$ID
hop.e <- hce_neg$ID
hop.f <- hcl_neg$ID

#GSE126108 Positive enhancers
hop.pos.enhancers <- list("DESeq2" = hop.a, "edgeR" = hop.b, "limmaVoom" = hop.c)
hpc.venn <- venn(hop.pos.enhancers, zcolor = "style")

hpc.int <- attributes(hpc.venn)
hpc.int <- hpc.int[["intersections"]][[7]]

#GSE126108 Negative enhancers
hop.neg.enhancers <- list("DESeq2" = hop.d, "edgeR" = hop.e, "limmaVoom" = hop.f)
hnc.venn <- venn(hop.neg.enhancers, zcolor = "style")

hnc.int <- attributes(hnc.venn)
hnc.int <- hnc.int[["intersections"]][[5]]

#####################################
#GSE127790 KO ES Cells Venn Diagrams #
#####################################

sh.a <- scd_pos$ID
sh.b <- sce_pos$ID
sh.c <- scl_pos$ID

sh.d <- scd_neg$ID
sh.e <- sce_neg$ID
sh.f <- scl_neg$ID


#GSE127790 Positive enhancers
sh.pos.enhancers <- list("DESeq2" = sh.a, "edgeR" = sh.b, "limmaVoom" = sh.c)
spc.venn <- venn(sh.pos.enhancers, zcolor = "style")

spc.int <- attributes(spc.venn)
spc.int <- spc.int[["intersections"]][[6]]

#GSE127790 Negative enhancers
sh.neg.enhancers <- list("DESeq2" = sh.d, "edgeR" = sh.e, "limmaVoom" = sh.f)
snc.venn <- venn(sh.neg.enhancers, zcolor = "style")

snc.int <- attributes(snc.venn)
snc.int <- snc.int[["intersections"]][[4]]


#Positive
#Venn diagram for GSE165689 and GSE126108 enhancer
nh.enhancer <- list("GSE165689 enhancer" = npc.int, "GSE126108 enhancer" = hpc.int)
nh.enhancer.venn <- venn(nh.enhancer, zcolor = "style")
nh.enhancer.int <- attributes(nh.enhancer.venn)
nh.enhancer.int <- nh.enhancer.int[["intersections"]][[3]]

#Venn diagram for GSE165689 and GSE127790 enhancer
ns.enhancer <- list("GSE165689 enhancer" = npc.int, "GSE127790 enhancer" = spc.int)
ns.enhancer.venn <- venn(ns.enhancer, zcolor = "style")
ns.enhancer.int <- attributes(ns.enhancer.venn)
ns.enhancer.int <- ns.enhancer.int[["intersections"]][[3]]

#Venn diagram for GSE165689, GSE126108, and GSE127790 enhancer
nhs.enhancer <- list("GSE165689 enhancer" = npc.int, "GSE127790 enhancer" = spc.int, "GSE126108 enhancer" = hpc.int)
nhs.enhancer.venn <- venn(nhs.enhancer, zcolor = "style", ilcs = 1)
nhs.enhancer.int <- attributes(nhs.enhancer.venn)
nhs.enhancer.int <- nhs.enhancer.int[["intersections"]][[6]]

#Negative
#Venn diagram for GSE165689 and GSE126108 enhancer
nh.enhancer <- list("GSE165689 enhancer" = nnc.int, "GSE126108 enhancer" = hnc.int)
nh.enhancer.venn <- venn(nh.enhancer, zcolor = "style")
nh.enhancer.int <- attributes(nh.enhancer.venn)
nh.enhancer.int <- nh.enhancer.int[["intersections"]][[]]

#Venn diagram for GSE165689 and GSE127790 enhancer
ns.enhancer <- list("GSE165689 enhancer" = nnc.int, "GSE127790 enhancer" = snc.int)
ns.enhancer.venn <- venn(ns.enhancer, zcolor = "style")
ns.enhancer.int <- attributes(ns.enhancer.venn)
ns.enhancer.int <- ns.enhancer.int[["intersections"]][[3]]

#Venn diagram for GSE165689, GSE126108, and GSE127790 enhancer
nhs.enhancer <- list("GSE165689 enhancer" = nnc.int, "GSE127790 enhancer" = snc.int, "GSE126108 enhancer" = hnc.int)
nhs.enhancer.venn <- venn(nhs.enhancer, zcolor = "style", ilcs = 1)
nhs.enhancer.int <- attributes(nhs.enhancer.venn)
nhs.enhancer.int <- nhs.enhancer.int[["intersections"]][[]]


#ocr
#########################################################################################################################################################
ncd <- read.delim("GSE165689/ZcKOvsWT_ocr_deseq2.txt")
nce <- read.delim("GSE165689/ZcKOvsWT_ocr_edgeR.txt")
ncl <- read.delim("GSE165689/ZcKOvsWT_ocr_limmaVoom.txt")

hcd <- read.delim("GSE126108 /ocr_deseq2.txt")
hce <- read.delim("GSE126108 /ocr_edgeR.txt")
hcl <- read.delim("GSE126108 /ocr_limmaVoom.txt")

scd <- read.delim("GSE127790/ocr_deseq2.txt")
sce <- read.delim("GSE127790/ocr_edgeR.txt")
scl <- read.delim("GSE127790/ocr_limmaVoom.txt")

ncd <- rownames_to_column(ncd, "ID")
nce <- rownames_to_column(nce, "ID")
ncl <- rownames_to_column(ncl, "ID")

hcd <- rownames_to_column(hcd, "ID")
hce <- rownames_to_column(hce, "ID")
hcl <- rownames_to_column(hcl, "ID")

scd <- rownames_to_column(scd, "ID")
sce <- rownames_to_column(sce, "ID")
scl <- rownames_to_column(scl, "ID")



#################################################
# Filter sets based on log2FC > 0.58 or < -0.58 #
#################################################
#GSE165689 Zcchc8 KO SIMS Cells

ncd_pos <- filter(ncd, ncd$log2FoldChange > 0.58 & ncd$padj <= 0.01)
ncd_neg <- filter(ncd, ncd$log2FoldChange < -0.58 & ncd$padj <= 0.01)

nce_pos <- filter(nce, nce$logFC > 0.58 & nce$FDR <= 0.01)
nce_neg <- filter(nce, nce$logFC < -0.58 & nce$FDR <= 0.01)

ncl_pos <- filter(ncl, ncl$logFC > 0.58 & ncl$adj.P.Val <= 0.01)
ncl_neg <- filter(ncl, ncl$logFC < -0.58 & ncl$adj.P.Val <= 0.01)


#GSE126108 KO mouse
hcd_pos <- filter(hcd, hcd$log2FoldChange > 0.58 & hcd$padj <= 0.01)
hcd_neg <- filter(hcd, hcd$log2FoldChange < -0.58 & hcd$padj <= 0.01)

hce_pos <- filter(hce, hce$logFC > 0.58 & hce$FDR <= 0.01)
hce_neg <- filter(hce, hce$logFC < -0.58 & hce$FDR <= 0.01)

hcl_pos <- filter(hcl, hcl$logFC > 0.58 & hcl$adj.P.Val <= 0.01)
hcl_neg <- filter(hcl, hcl$logFC < -0.58 & hcl$adj.P.Val <= 0.01)


#GSE127790 KO ES cells
scd_pos <- filter(scd, scd$log2FoldChange > 0.58 & scd$padj <= 0.01)
scd_neg <- filter(scd, scd$log2FoldChange < -0.58 & scd$padj <= 0.01)

sce_pos <- filter(sce, sce$logFC > 0.58 & sce$FDR <= 0.01)
sce_neg <- filter(sce, sce$logFC < -0.58 & sce$FDR <= 0.01)

scl_pos <- filter(scl, scl$logFC > 0.58 & scl$adj.P.Val <= 0.01)
scl_neg <- filter(scl, scl$logFC < -0.58 & scl$adj.P.Val <= 0.01)


################################
# GSE165689 SIMS Cells Venn Diagrams #
################################

a <- ncd_pos$ID
b <- nce_pos$ID
c <- ncl_pos$ID

d <- ncd_neg$ID
e <- nce_neg$ID
f <- ncl_neg$ID


#GSE165689 Positive ocr
GSE165689.pos.ocr <- list("DESeq2" = a, "edgeR" = b, "Limma-Voom" = c)
npc.venn <- venn(GSE165689.pos.ocr, zcolor = "style", ilcs = 1, sncs = 1.2)

npc.int <- attributes(npc.venn)
npc.int <- npc.int[["intersections"]][[6]]


#GSE165689 Negative ocr
GSE165689.neg.ocr <- list("DESeq2" = d, "edgeR" = e, "limmaVoom" = f)
nnc.venn <- venn(GSE165689.neg.ocr, zcolor = "style")

nnc.int <- attributes(nnc.venn)
nnc.int <- nnc.int[["intersections"]][[6]]


#################################
# GSE126108 KO MouseVenn Diagrams #
#################################

hop.a <- hcd_pos$ID
hop.b <- hce_pos$ID
hop.c <- hcl_pos$ID

hop.d <- hcd_neg$ID
hop.e <- hce_neg$ID
hop.f <- hcl_neg$ID

#GSE126108 Positive ocr
hop.pos.ocr <- list("DESeq2" = hop.a, "edgeR" = hop.b, "limmaVoom" = hop.c)
hpc.venn <- venn(hop.pos.ocr, zcolor = "style")

hpc.int <- attributes(hpc.venn)
hpc.int <- hpc.int[["intersections"]][[7]]

#GSE126108 Negative ocr
hop.neg.ocr <- list("DESeq2" = hop.d, "edgeR" = hop.e, "limmaVoom" = hop.f)
hnc.venn <- venn(hop.neg.ocr, zcolor = "style")

hnc.int <- attributes(hnc.venn)
hnc.int <- hnc.int[["intersections"]][[5]]

#####################################
#GSE127790 KO ES Cells Venn Diagrams #
#####################################

sh.a <- scd_pos$ID
sh.b <- sce_pos$ID
sh.c <- scl_pos$ID

sh.d <- scd_neg$ID
sh.e <- sce_neg$ID
sh.f <- scl_neg$ID


#GSE127790 Positive ocr
sh.pos.ocr <- list("DESeq2" = sh.a, "edgeR" = sh.b, "limmaVoom" = sh.c)
spc.venn <- venn(sh.pos.ocr, zcolor = "style")

spc.int <- attributes(spc.venn)
spc.int <- spc.int[["intersections"]][[7]]

#GSE127790 Negative ocr
sh.neg.ocr <- list("DESeq2" = sh.d, "edgeR" = sh.e, "limmaVoom" = sh.f)
snc.venn <- venn(sh.neg.ocr, zcolor = "style")

snc.int <- attributes(snc.venn)
snc.int <- snc.int[["intersections"]][[4]]


#Positive
#Venn diagram for GSE165689 and GSE126108 ocr
nh.ocr <- list("GSE165689 ocr" = npc.int, "GSE126108 ocr" = hpc.int)
nh.ocr.venn <- venn(nh.ocr, zcolor = "style")
nh.ocr.int <- attributes(nh.ocr.venn)
nh.ocr.int <- nh.ocr.int[["intersections"]][[3]]

#Venn diagram for GSE165689 and GSE127790 ocr
ns.ocr <- list("GSE165689 ocr" = npc.int, "GSE127790 ocr" = spc.int)
ns.ocr.venn <- venn(ns.ocr, zcolor = "style")
ns.ocr.int <- attributes(ns.ocr.venn)
ns.ocr.int <- ns.ocr.int[["intersections"]][[3]]

#Venn diagram for GSE165689, GSE126108, and GSE127790 ocr
nhs.ocr <- list("GSE165689 ocr" = npc.int, "GSE127790 ocr" = spc.int, "GSE126108 ocr" = hpc.int)
nhs.ocr.venn <- venn(nhs.ocr, zcolor = "style", ilcs = 1)
nhs.ocr.int <- attributes(nhs.ocr.venn)
nhs.ocr.int <- nhs.ocr.int[["intersections"]][[6]]

#Negative
#Venn diagram for GSE165689 and GSE126108 ocr
nh.ocr <- list("GSE165689 ocr" = nnc.int, "GSE126108 ocr" = hnc.int)
nh.ocr.venn <- venn(nh.ocr, zcolor = "style")
nh.ocr.int <- attributes(nh.ocr.venn)
nh.ocr.int <- nh.ocr.int[["intersections"]][[]]

#Venn diagram for GSE165689 and GSE127790 ocr
ns.ocr <- list("GSE165689 ocr" = nnc.int, "GSE127790 ocr" = snc.int)
ns.ocr.venn <- venn(ns.ocr, zcolor = "style")
ns.ocr.int <- attributes(ns.ocr.venn)
ns.ocr.int <- ns.ocr.int[["intersections"]][[3]]

#Venn diagram for GSE165689, GSE126108, and GSE127790 ocr
nhs.ocr <- list("GSE165689 ocr" = nnc.int, "GSE127790 ocr" = snc.int, "GSE126108 ocr" = hnc.int)
nhs.ocr.venn <- venn(nhs.ocr, zcolor = "style", ilcs = 1)
nhs.ocr.int <- attributes(nhs.ocr.venn)
nhs.ocr.int <- nhs.ocr.int[["intersections"]][[]]


#promoters
#########################################################################################################################################################
ncd <- read.delim("GSE165689/ZcKOvsWT_promoters_deseq2.txt")
nce <- read.delim("GSE165689/ZcKOvsWT_promoters_edgeR.txt")
ncl <- read.delim("GSE165689/ZcKOvsWT_promoters_limmaVoom.txt")

hcd <- read.delim("GSE126108 /promoters_deseq2.txt")
hce <- read.delim("GSE126108 /promoters_edgeR.txt")
hcl <- read.delim("GSE126108 /promoters_limmaVoom.txt")

scd <- read.delim("GSE127790/promoters/promoters_deseq2.txt")
sce <- read.delim("GSE127790/promoters/promoters_edgeR.txt")
scl <- read.delim("GSE127790/promoters/promoters_limmaVoom.txt")

ncd <- rownames_to_column(ncd, "ID")
nce <- rownames_to_column(nce, "ID")
ncl <- rownames_to_column(ncl, "ID")

hcd <- rownames_to_column(hcd, "ID")
hce <- rownames_to_column(hce, "ID")
hcl <- rownames_to_column(hcl, "ID")

scd <- rownames_to_column(scd, "ID")
sce <- rownames_to_column(sce, "ID")
scl <- rownames_to_column(scl, "ID")



#################################################
# Filter sets based on log2FC > 0.58 or < -0.58 #
#################################################
#GSE165689 Zcchc8 KO SIMS Cells

ncd_pos <- filter(ncd, ncd$log2FoldChange > 0.58 & ncd$padj <= 0.01)
ncd_neg <- filter(ncd, ncd$log2FoldChange < -0.58 & ncd$padj <= 0.01)

nce_pos <- filter(nce, nce$logFC > 0.58 & nce$FDR <= 0.01)
nce_neg <- filter(nce, nce$logFC < -0.58 & nce$FDR <= 0.01)

ncl_pos <- filter(ncl, ncl$logFC > 0.58 & ncl$adj.P.Val <= 0.01)
ncl_neg <- filter(ncl, ncl$logFC < -0.58 & ncl$adj.P.Val <= 0.01)


#GSE126108 KO mouse
hcd_pos <- filter(hcd, hcd$log2FoldChange > 0.58 & hcd$padj <= 0.01)
hcd_neg <- filter(hcd, hcd$log2FoldChange < -0.58 & hcd$padj <= 0.01)

hce_pos <- filter(hce, hce$logFC > 0.58 & hce$FDR <= 0.01)
hce_neg <- filter(hce, hce$logFC < -0.58 & hce$FDR <= 0.01)

hcl_pos <- filter(hcl, hcl$logFC > 0.58 & hcl$adj.P.Val <= 0.01)
hcl_neg <- filter(hcl, hcl$logFC < -0.58 & hcl$adj.P.Val <= 0.01)


#GSE127790 KO ES cells
scd_pos <- filter(scd, scd$log2FoldChange > 0.58 & scd$padj <= 0.01)
scd_neg <- filter(scd, scd$log2FoldChange < -0.58 & scd$padj <= 0.01)

sce_pos <- filter(sce, sce$logFC > 0.58 & sce$FDR <= 0.01)
sce_neg <- filter(sce, sce$logFC < -0.58 & sce$FDR <= 0.01)

scl_pos <- filter(scl, scl$logFC > 0.58 & scl$adj.P.Val <= 0.01)
scl_neg <- filter(scl, scl$logFC < -0.58 & scl$adj.P.Val <= 0.01)


################################
# GSE165689 SIMS Cells Venn Diagrams #
################################

a <- ncd_pos$ID
b <- nce_pos$ID
c <- ncl_pos$ID

d <- ncd_neg$ID
e <- nce_neg$ID
f <- ncl_neg$ID


#GSE165689 Positive promoters
GSE165689.pos.promoters <- list("DESeq2" = a, "edgeR" = b, "Limma-Voom" = c)
npc.venn <- venn(GSE165689.pos.promoters, zcolor = "style", ilcs = 1, sncs = 1.2)

npc.int <- attributes(npc.venn)
npc.int <- npc.int[["intersections"]][[6]]


#GSE165689 Negative promoters
GSE165689.neg.promoters <- list("DESeq2" = d, "edgeR" = e, "limmaVoom" = f)
nnc.venn <- venn(GSE165689.neg.promoters, zcolor = "style")

nnc.int <- attributes(nnc.venn)
nnc.int <- nnc.int[["intersections"]][[5]]


#################################
# GSE126108 KO MouseVenn Diagrams #
#################################

hop.a <- hcd_pos$ID
hop.b <- hce_pos$ID
hop.c <- hcl_pos$ID

hop.d <- hcd_neg$ID
hop.e <- hce_neg$ID
hop.f <- hcl_neg$ID

#GSE126108 Positive promoters
hop.pos.promoters <- list("DESeq2" = hop.a, "edgeR" = hop.b, "limmaVoom" = hop.c)
hpc.venn <- venn(hop.pos.promoters, zcolor = "style")

hpc.int <- attributes(hpc.venn)
hpc.int <- hpc.int[["intersections"]][[5]]

#GSE126108 Negative promoters
hop.neg.promoters <- list("DESeq2" = hop.d, "edgeR" = hop.e, "limmaVoom" = hop.f)
hnc.venn <- venn(hop.neg.promoters, zcolor = "style")

hnc.int <- attributes(hnc.venn)
hnc.int <- hnc.int[["intersections"]][[6]]

#####################################
#GSE127790 KO ES Cells Venn Diagrams #
#####################################

sh.a <- scd_pos$ID
sh.b <- sce_pos$ID
sh.c <- scl_pos$ID

sh.d <- scd_neg$ID
sh.e <- sce_neg$ID
sh.f <- scl_neg$ID


#GSE127790 Positive promoters
sh.pos.promoters <- list("DESeq2" = sh.a, "edgeR" = sh.b, "limmaVoom" = sh.c)
spc.venn <- venn(sh.pos.promoters, zcolor = "style")

spc.int <- attributes(spc.venn)
spc.int <- spc.int[["intersections"]][[6]]

#GSE127790 Negative promoters
sh.neg.promoters <- list("DESeq2" = sh.d, "edgeR" = sh.e, "limmaVoom" = sh.f)
snc.venn <- venn(sh.neg.promoters, zcolor = "style")

snc.int <- attributes(snc.venn)
snc.int <- snc.int[["intersections"]][[5]]


#Positive
#Venn diagram for GSE165689 and GSE126108 promoters
nh.promoters <- list("GSE165689 promoters" = npc.int, "GSE126108 promoters" = hpc.int)
nh.promoters.venn <- venn(nh.promoters, zcolor = "style")
nh.promoters.int <- attributes(nh.promoters.venn)
nh.promoters.int <- nh.promoters.int[["intersections"]][[3]]

#Venn diagram for GSE165689 and GSE127790 promoters
ns.promoters <- list("GSE165689 promoters" = npc.int, "GSE127790 promoters" = spc.int)
ns.promoters.venn <- venn(ns.promoters, zcolor = "style")
ns.promoters.int <- attributes(ns.promoters.venn)
ns.promoters.int <- ns.promoters.int[["intersections"]][[3]]

#Venn diagram for GSE165689, GSE126108, and GSE127790 promoters
nhs.promoters <- list("GSE165689 promoters" = npc.int, "GSE127790 promoters" = spc.int, "GSE126108 promoters" = hpc.int)
nhs.promoters.venn <- venn(nhs.promoters, zcolor = "style", ilcs = 1)
nhs.promoters.int <- attributes(nhs.promoters.venn)
nhs.promoters.int <- nhs.promoters.int[["intersections"]][[7]]

#Negative
#Venn diagram for GSE165689 and GSE126108 promoters
nh.promoters <- list("GSE165689 promoters" = nnc.int, "GSE126108 promoters" = hnc.int)
nh.promoters.venn <- venn(nh.promoters, zcolor = "style")
nh.promoters.int <- attributes(nh.promoters.venn)
nh.promoters.int <- nh.promoters.int[["intersections"]][[]]

#Venn diagram for GSE165689 and GSE127790 promoters
ns.promoters <- list("GSE165689 promoters" = nnc.int, "GSE127790 promoters" = snc.int)
ns.promoters.venn <- venn(ns.promoters, zcolor = "style")
ns.promoters.int <- attributes(ns.promoters.venn)
ns.promoters.int <- ns.promoters.int[["intersections"]][[3]]

#Venn diagram for GSE165689, GSE126108, and GSE127790 promoters
nhs.promoters <- list("GSE165689 promoters" = nnc.int, "GSE127790 promoters" = snc.int, "GSE126108 promoters" = hnc.int)
nhs.promoters.venn <- venn(nhs.promoters, zcolor = "style", ilcs = 1)
nhs.promoters.int <- attributes(nhs.promoters.venn)
nhs.promoters.int <- nhs.promoters.int[["intersections"]][[]]


#pfr
#########################################################################################################################################################
ncd <- read.delim("GSE165689/ZcKOvsWT_pfr_deseq2.txt")
nce <- read.delim("GSE165689/ZcKOvsWT_pfr_edgeR.txt")
ncl <- read.delim("GSE165689/ZcKOvsWT_pfr_limmaVoom.txt")

hcd <- read.delim("GSE126108 /pfr_deseq2.txt")
hce <- read.delim("GSE126108 /pfr_edgeR.txt")
hcl <- read.delim("GSE126108 /pfr_limmaVoom.txt")

scd <- read.delim("GSE127790/pfr/pfr_deseq2.txt")
sce <- read.delim("GSE127790/pfr/pfr_edgeR.txt")
scl <- read.delim("GSE127790/pfr/pfr_limmaVoom.txt")

ncd <- rownames_to_column(ncd, "ID")
nce <- rownames_to_column(nce, "ID")
ncl <- rownames_to_column(ncl, "ID")

hcd <- rownames_to_column(hcd, "ID")
hce <- rownames_to_column(hce, "ID")
hcl <- rownames_to_column(hcl, "ID")

scd <- rownames_to_column(scd, "ID")
sce <- rownames_to_column(sce, "ID")
scl <- rownames_to_column(scl, "ID")



#################################################
# Filter sets based on log2FC > 0.58 or < -0.58 #
#################################################
#GSE165689 Zcchc8 KO SIMS Cells

ncd_pos <- filter(ncd, ncd$log2FoldChange > 0.58 & ncd$padj <= 0.01)
ncd_neg <- filter(ncd, ncd$log2FoldChange < -0.58 & ncd$padj <= 0.01)

nce_pos <- filter(nce, nce$logFC > 0.58 & nce$FDR <= 0.01)
nce_neg <- filter(nce, nce$logFC < -0.58 & nce$FDR <= 0.01)

ncl_pos <- filter(ncl, ncl$logFC > 0.58 & ncl$adj.P.Val <= 0.01)
ncl_neg <- filter(ncl, ncl$logFC < -0.58 & ncl$adj.P.Val <= 0.01)


#GSE126108 KO mouse
hcd_pos <- filter(hcd, hcd$log2FoldChange > 0.58 & hcd$padj <= 0.01)
hcd_neg <- filter(hcd, hcd$log2FoldChange < -0.58 & hcd$padj <= 0.01)

hce_pos <- filter(hce, hce$logFC > 0.58 & hce$FDR <= 0.01)
hce_neg <- filter(hce, hce$logFC < -0.58 & hce$FDR <= 0.01)

hcl_pos <- filter(hcl, hcl$logFC > 0.58 & hcl$adj.P.Val <= 0.01)
hcl_neg <- filter(hcl, hcl$logFC < -0.58 & hcl$adj.P.Val <= 0.01)


#GSE127790 KO ES cells
scd_pos <- filter(scd, scd$log2FoldChange > 0.58 & scd$padj <= 0.01)
scd_neg <- filter(scd, scd$log2FoldChange < -0.58 & scd$padj <= 0.01)

sce_pos <- filter(sce, sce$logFC > 0.58 & sce$FDR <= 0.01)
sce_neg <- filter(sce, sce$logFC < -0.58 & sce$FDR <= 0.01)

scl_pos <- filter(scl, scl$logFC > 0.58 & scl$adj.P.Val <= 0.01)
scl_neg <- filter(scl, scl$logFC < -0.58 & scl$adj.P.Val <= 0.01)


################################
# GSE165689 SIMS Cells Venn Diagrams #
################################

a <- ncd_pos$ID
b <- nce_pos$ID
c <- ncl_pos$ID

d <- ncd_neg$ID
e <- nce_neg$ID
f <- ncl_neg$ID


#GSE165689 Positive pfr
GSE165689.pos.pfr <- list("DESeq2" = a, "edgeR" = b, "Limma-Voom" = c)
npc.venn <- venn(GSE165689.pos.pfr, zcolor = "style", ilcs = 1, sncs = 1.2)

npc.int <- attributes(npc.venn)
npc.int <- npc.int[["intersections"]][[7]]


#GSE165689 Negative pfr
GSE165689.neg.pfr <- list("DESeq2" = d, "edgeR" = e, "limmaVoom" = f)
nnc.venn <- venn(GSE165689.neg.pfr, zcolor = "style")

nnc.int <- attributes(nnc.venn)
nnc.int <- nnc.int[["intersections"]][[5]]


#################################
# GSE126108 KO MouseVenn Diagrams #
#################################

hop.a <- hcd_pos$ID
hop.b <- hce_pos$ID
hop.c <- hcl_pos$ID

hop.d <- hcd_neg$ID
hop.e <- hce_neg$ID
hop.f <- hcl_neg$ID

#GSE126108 Positive pfr
hop.pos.pfr <- list("DESeq2" = hop.a, "edgeR" = hop.b, "limmaVoom" = hop.c)
hpc.venn <- venn(hop.pos.pfr, zcolor = "style")

hpc.int <- attributes(hpc.venn)
hpc.int <- hpc.int[["intersections"]][[7]]

#GSE126108 Negative pfr
hop.neg.pfr <- list("DESeq2" = hop.d, "edgeR" = hop.e, "limmaVoom" = hop.f)
hnc.venn <- venn(hop.neg.pfr, zcolor = "style")

hnc.int <- attributes(hnc.venn)
hnc.int <- hnc.int[["intersections"]][[6]]

#####################################
#GSE127790 KO ES Cells Venn Diagrams #
#####################################

sh.a <- scd_pos$ID
sh.b <- sce_pos$ID
sh.c <- scl_pos$ID

sh.d <- scd_neg$ID
sh.e <- sce_neg$ID
sh.f <- scl_neg$ID


#GSE127790 Positive pfr
sh.pos.pfr <- list("DESeq2" = sh.a, "edgeR" = sh.b, "limmaVoom" = sh.c)
spc.venn <- venn(sh.pos.pfr, zcolor = "style")

spc.int <- attributes(spc.venn)
spc.int <- spc.int[["intersections"]][[6]]

#GSE127790 Negative pfr
sh.neg.pfr <- list("DESeq2" = sh.d, "edgeR" = sh.e, "limmaVoom" = sh.f)
snc.venn <- venn(sh.neg.pfr, zcolor = "style")

snc.int <- attributes(snc.venn)
snc.int <- snc.int[["intersections"]][[5]]


#Positive
#Venn diagram for GSE165689 and GSE126108 pfr
nh.pfr <- list("GSE165689 pfr" = npc.int, "GSE126108 pfr" = hpc.int)
nh.pfr.venn <- venn(nh.pfr, zcolor = "style")
nh.pfr.int <- attributes(nh.pfr.venn)
nh.pfr.int <- nh.pfr.int[["intersections"]][[3]]

#Venn diagram for GSE165689 and GSE127790 pfr
ns.pfr <- list("GSE165689 pfr" = npc.int, "GSE127790 pfr" = spc.int)
ns.pfr.venn <- venn(ns.pfr, zcolor = "style")
ns.pfr.int <- attributes(ns.pfr.venn)
ns.pfr.int <- ns.pfr.int[["intersections"]][[3]]

#Venn diagram for GSE165689, GSE126108, and GSE127790 pfr
nhs.pfr <- list("GSE165689 pfr" = npc.int, "GSE127790 pfr" = spc.int, "GSE126108 pfr" = hpc.int)
nhs.pfr.venn <- venn(nhs.pfr, zcolor = "style", ilcs = 1)
nhs.pfr.int <- attributes(nhs.pfr.venn)
nhs.pfr.int <- nhs.pfr.int[["intersections"]][[7]]

#Negative
#Venn diagram for GSE165689 and GSE126108 pfr
nh.pfr <- list("GSE165689 pfr" = nnc.int, "GSE126108 pfr" = hnc.int)
nh.pfr.venn <- venn(nh.pfr, zcolor = "style")
nh.pfr.int <- attributes(nh.pfr.venn)
nh.pfr.int <- nh.pfr.int[["intersections"]][[]]

#Venn diagram for GSE165689 and GSE127790 pfr
ns.pfr <- list("GSE165689 pfr" = nnc.int, "GSE127790 pfr" = snc.int)
ns.pfr.venn <- venn(ns.pfr, zcolor = "style")
ns.pfr.int <- attributes(ns.pfr.venn)
ns.pfr.int <- ns.pfr.int[["intersections"]][[3]]

#Venn diagram for GSE165689, GSE126108, and GSE127790 pfr
nhs.pfr <- list("GSE165689 pfr" = nnc.int, "GSE127790 pfr" = snc.int, "GSE126108 pfr" = hnc.int)
nhs.pfr.venn <- venn(nhs.pfr, zcolor = "style", ilcs = 1)
nhs.pfr.int <- attributes(nhs.pfr.venn)
nhs.pfr.int <- nhs.pfr.int[["intersections"]][[]]


#tfbs
#########################################################################################################################################################
ncd <- read.delim("GSE165689/ZcKOvsWT_tfbs_deseq2.txt")
nce <- read.delim("GSE165689/ZcKOvsWT_tfbs_edgeR.txt")
ncl <- read.delim("GSE165689/ZcKOvsWT_tfbs_limmaVoom.txt")

hcd <- read.delim("GSE126108 /tfbs_deseq2.txt")
hce <- read.delim("GSE126108 /tfbs_edgeR.txt")
hcl <- read.delim("GSE126108 /tfbs_limmaVoom.txt")

scd <- read.delim("GSE127790/tfbs_deseq2.txt")
sce <- read.delim("GSE127790/tfbs_edgeR.txt")
scl <- read.delim("GSE127790/tfbs_limmaVoom.txt")

ncd <- rownames_to_column(ncd, "ID")
nce <- rownames_to_column(nce, "ID")
ncl <- rownames_to_column(ncl, "ID")

hcd <- rownames_to_column(hcd, "ID")
hce <- rownames_to_column(hce, "ID")
hcl <- rownames_to_column(hcl, "ID")

scd <- rownames_to_column(scd, "ID")
sce <- rownames_to_column(sce, "ID")
scl <- rownames_to_column(scl, "ID")



#################################################
# Filter sets based on log2FC > 0.58 or < -0.58 #
#################################################
#GSE165689 Zcchc8 KO SIMS Cells

ncd_pos <- filter(ncd, ncd$log2FoldChange > 0.58 & ncd$padj <= 0.01)
ncd_neg <- filter(ncd, ncd$log2FoldChange < -0.58 & ncd$padj <= 0.01)

nce_pos <- filter(nce, nce$logFC > 0.58 & nce$FDR <= 0.01)
nce_neg <- filter(nce, nce$logFC < -0.58 & nce$FDR <= 0.01)

ncl_pos <- filter(ncl, ncl$logFC > 0.58 & ncl$adj.P.Val <= 0.01)
ncl_neg <- filter(ncl, ncl$logFC < -0.58 & ncl$adj.P.Val <= 0.01)


#GSE126108 KO mouse
hcd_pos <- filter(hcd, hcd$log2FoldChange > 0.58 & hcd$padj <= 0.01)
hcd_neg <- filter(hcd, hcd$log2FoldChange < -0.58 & hcd$padj <= 0.01)

hce_pos <- filter(hce, hce$logFC > 0.58 & hce$FDR <= 0.01)
hce_neg <- filter(hce, hce$logFC < -0.58 & hce$FDR <= 0.01)

hcl_pos <- filter(hcl, hcl$logFC > 0.58 & hcl$adj.P.Val <= 0.01)
hcl_neg <- filter(hcl, hcl$logFC < -0.58 & hcl$adj.P.Val <= 0.01)


#GSE127790 KO ES cells
scd_pos <- filter(scd, scd$log2FoldChange > 0.58 & scd$padj <= 0.01)
scd_neg <- filter(scd, scd$log2FoldChange < -0.58 & scd$padj <= 0.01)

sce_pos <- filter(sce, sce$logFC > 0.58 & sce$FDR <= 0.01)
sce_neg <- filter(sce, sce$logFC < -0.58 & sce$FDR <= 0.01)

scl_pos <- filter(scl, scl$logFC > 0.58 & scl$adj.P.Val <= 0.01)
scl_neg <- filter(scl, scl$logFC < -0.58 & scl$adj.P.Val <= 0.01)


################################
# GSE165689 SIMS Cells Venn Diagrams #
################################

a <- ncd_pos$ID
b <- nce_pos$ID
c <- ncl_pos$ID

d <- ncd_neg$ID
e <- nce_neg$ID
f <- ncl_neg$ID


#GSE165689 Positive tfbs
GSE165689.pos.tfbs <- list("DESeq2" = a, "edgeR" = b, "Limma-Voom" = c)
npc.venn <- venn(GSE165689.pos.tfbs, zcolor = "style", ilcs = 1, sncs = 1.2)

npc.int <- attributes(npc.venn)
npc.int <- npc.int[["intersections"]][[5]]


#GSE165689 Negative tfbs
GSE165689.neg.tfbs <- list("DESeq2" = d, "edgeR" = e, "limmaVoom" = f)
nnc.venn <- venn(GSE165689.neg.tfbs, zcolor = "style")

nnc.int <- attributes(nnc.venn)
nnc.int <- nnc.int[["intersections"]][[6]]


#################################
# GSE126108 KO MouseVenn Diagrams #
#################################

hop.a <- hcd_pos$ID
hop.b <- hce_pos$ID
hop.c <- hcl_pos$ID

hop.d <- hcd_neg$ID
hop.e <- hce_neg$ID
hop.f <- hcl_neg$ID

#GSE126108 Positive tfbs
hop.pos.tfbs <- list("DESeq2" = hop.a, "edgeR" = hop.b, "limmaVoom" = hop.c)
hpc.venn <- venn(hop.pos.tfbs, zcolor = "style")

hpc.int <- attributes(hpc.venn)
hpc.int <- hpc.int[["intersections"]][[5]]

#GSE126108 Negative tfbs
hop.neg.tfbs <- list("DESeq2" = hop.d, "edgeR" = hop.e, "limmaVoom" = hop.f)
hnc.venn <- venn(hop.neg.tfbs, zcolor = "style")

hnc.int <- attributes(hnc.venn)
hnc.int <- hnc.int[["intersections"]][[5]]

#####################################
#GSE127790 KO ES Cells Venn Diagrams #
#####################################

sh.a <- scd_pos$ID
sh.b <- sce_pos$ID
sh.c <- scl_pos$ID

sh.d <- scd_neg$ID
sh.e <- sce_neg$ID
sh.f <- scl_neg$ID


#GSE127790 Positive tfbs
sh.pos.tfbs <- list("DESeq2" = sh.a, "edgeR" = sh.b, "limmaVoom" = sh.c)
spc.venn <- venn(sh.pos.tfbs, zcolor = "style")

spc.int <- attributes(spc.venn)
spc.int <- spc.int[["intersections"]][[4]]

#GSE127790 Negative tfbs
sh.neg.tfbs <- list("DESeq2" = sh.d, "edgeR" = sh.e, "limmaVoom" = sh.f)
snc.venn <- venn(sh.neg.tfbs, zcolor = "style")

snc.int <- attributes(snc.venn)
snc.int <- snc.int[["intersections"]][[5]]


#Positive
#Venn diagram for GSE165689 and GSE126108 tfbs
nh.tfbs <- list("GSE165689 tfbs" = npc.int, "GSE126108 tfbs" = hpc.int)
nh.tfbs.venn <- venn(nh.tfbs, zcolor = "style")
nh.tfbs.int <- attributes(nh.tfbs.venn)
nh.tfbs.int <- nh.tfbs.int[["intersections"]][[3]]

#Venn diagram for GSE165689 and GSE127790 tfbs
ns.tfbs <- list("GSE165689 tfbs" = npc.int, "GSE127790 tfbs" = spc.int)
ns.tfbs.venn <- venn(ns.tfbs, zcolor = "style")
ns.tfbs.int <- attributes(ns.tfbs.venn)
ns.tfbs.int <- ns.tfbs.int[["intersections"]][[3]]

#Venn diagram for GSE165689, GSE126108, and GSE127790 tfbs
nhs.tfbs <- list("GSE165689 tfbs" = npc.int, "GSE127790 tfbs" = spc.int, "GSE126108 tfbs" = hpc.int)
nhs.tfbs.venn <- venn(nhs.tfbs, zcolor = "style", ilcs = 1)
nhs.tfbs.int <- attributes(nhs.tfbs.venn)
nhs.tfbs.int <- nhs.tfbs.int[["intersections"]][[7]]

#Negative
#Venn diagram for GSE165689 and GSE126108 tfbs
nh.tfbs <- list("GSE165689 tfbs" = nnc.int, "GSE126108 tfbs" = hnc.int)
nh.tfbs.venn <- venn(nh.tfbs, zcolor = "style")
nh.tfbs.int <- attributes(nh.tfbs.venn)
nh.tfbs.int <- nh.tfbs.int[["intersections"]][[]]

#Venn diagram for GSE165689 and GSE127790 tfbs
ns.tfbs <- list("GSE165689 tfbs" = nnc.int, "GSE127790 tfbs" = snc.int)
ns.tfbs.venn <- venn(ns.tfbs, zcolor = "style")
ns.tfbs.int <- attributes(ns.tfbs.venn)
ns.tfbs.int <- ns.tfbs.int[["intersections"]][[3]]

#Venn diagram for GSE165689, GSE126108, and GSE127790 tfbs
nhs.tfbs <- list("GSE165689 tfbs" = nnc.int, "GSE127790 tfbs" = snc.int, "GSE126108 tfbs" = hnc.int)
nhs.tfbs.venn <- venn(nhs.tfbs, zcolor = "style", ilcs = 1)
nhs.tfbs.int <- attributes(nhs.tfbs.venn)
nhs.tfbs.int <- nhs.tfbs.int[["intersections"]][[]]